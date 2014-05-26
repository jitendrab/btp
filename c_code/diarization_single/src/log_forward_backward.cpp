/*
Written By: Jitendra Prasad Keer
BTech, CSE, IIT Mandi
 */

/******************************************************************************
One forward-backward pass through a minimum-duration HMM model with a
single Gaussian in each of the states. T: totalFeatures
*******************************************************************************/
double ESHMM::mdHMMLogForwardBackward(ESHMM *mdHMM, VECTOR_OF_F_VECTORS *features, double **post, int T, mat &gamma, rowvec &gamma1,
			     mat &sumxi){
  printf("forward backward algorithm calculation in progress...\n");
  int i = 0, j = 0 , k = 0;
  /* total no of states is much larger, instead of number of pdfs we have to extend 
   states by Min_DUR, therefore total states = Q * MD */
  int Q = mdHMM->hmmStates;
  int Qmd = Q * MIN_DUR;
  mat logalpha(T, Qmd); // forward probability matrix
  mat logbeta(T, Qmd); // backward probability matrix
  mat logg(T, Qmd);  // loggamma 
  mat m(Q, 1);
  mat logA(Qmd, Qmd);  /// transition matrix is already in logarithm
  mat new_logp(T, Qmd); // after replication for each substates
  mat logp_k(Q, T);  // we have single cluster only, probability of each feature corresponding to each cluster
  
  printf("Q: %d  Qmd: %d\n", Q, Qmd);
  for(i = 0; i < Qmd; i++){
    for(j = 0; j < Qmd; j++){
      logA(i, j) = mdHMM->trans[i]->array[j];
    }
  }
  
  // minimum duration viterbi hence modify B(posterior) prob matrix
  
  for(i = 0; i < Q; i++)
    for(j = 0; j < T; j++){
      logp_k(i, j) = post[i][j];
    }
  
  for(i = 0; i < Q; i++){
    m(i, 0) = 1;
    for(j = 0; j < T; j++)
      logp_k(i, j) = 0.0; // since  we have only one cluster so cluster probability and 
    // total probability is same. Hence subtracting cluster probability from total probability would make it zero.
  }
  
  // modifying logp matrix according to minimum duration
  for(i = 0; i < Q; i++){
    for(j = 0; j < T; j++){
      for(k = i*MIN_DUR; k < (i+1)*MIN_DUR; k++){
	new_logp(j, k) = post[i][j];
      }
    }
  }
  /* forward initialization */
  // for summing log probabilties, first sum probs and then take logarithm
  printf("forward initialization...\n\n");
  for(i = 0; i < Qmd; i++){
    logalpha(0, i) = mdHMM->prior->array[i] + new_logp(0, i) ;
  }
  ///print logalpha after initialization
  for(i = 0; i < Qmd; i++)
    printf("%lf ", logalpha(0, i));
  
  /* forward induction */
  printf("forward induction in progress...\n");
  int t = 0;  
  
  mpfr_t summation3;
  mpfr_init(summation3);
  mpfr_t var11, var21;
  mpfr_init(var11);
  mpfr_init(var21);
  mpfr_set_d(var11, 0.0, MPFR_RNDN);
  mpfr_set_d(var21, 0.0, MPFR_RNDN);
  mpfr_set_d(summation3, 0.0, MPFR_RNDN);
  
  for(t = 1; t < T; t++){
    //printf("%d ", t);
    for(j = 0; j < Qmd; j++){
      vec v1(Qmd), v2(Qmd);      vec v3(Qmd);  
      //first find logalpha vector
      for(i = 0; i < Qmd; i++)
	v1(i) = logalpha(t-1, i);
      // if(t < 20)
      // 	v1.print("v1:\n");
      
      // extract transition probability vector
      for(i = 0; i < Qmd; i++)
	v2(i) = logA(i, j);
      // if(t < 20)
      // 	v2.print("v2:\n");
      // Now sum both the vectors into one
      for(i = 0; i < Qmd; i++)
	v3(i) = v1(i) + v2(i);
            
      double *temp = (double *)calloc(Qmd, sizeof(double ));
      for(i = 0; i < Qmd; i++)
	temp[i] = v3(i);
      // if(t < 20)
      // 	v3.print("v3:\n");
      //printf("printed\n");
      // now sum over whole column vector      
      mpfr_set_d(summation3, 0.0, MPFR_RNDN);
      // take the exponentiation and summation in one loop 
      
      // getting double from mpfr variable
      /// double mpfr_get_d(mpfr_t op, mpfr_rnd_t rnd);      
      //mpfr_set_d(var1, 0.0, MPFR_RNDD);
      //mpfr_set_d(var2, 0.0, MPFR_RNDD);
      // now take the exponentiation
      for(i = 0; i < Qmd; i++){
	double elem = temp[i];
	mpfr_set_d(var21, elem, MPFR_RNDD);
	//mpfr_printf("var2: %lf\n", var21);
	mpfr_exp(var11, var21, MPFR_RNDD); ///take exp(v2) and store in v1 
	// take sum of all elements in total  
	mpfr_add(summation3, summation3, var11, MPFR_RNDD); // add summation and v1
      }
      
      // now take the logarithm of sum 
      mpfr_log(summation3, summation3, MPFR_RNDD);
      // now convert this sum to double
      double sum2 = mpfr_get_d(summation3, MPFR_RNDD);
      // now assign this double to logalpha
      
      // now add logp(t, j)      
      sum2 += new_logp(t, j);
      // if(t < 20)
      // 	printf("sum: %lf\n", sum2);
      logalpha(t, j) = sum2;
      /// clear mpfr variables      
    }
    if(t < 20){
      printf("logalpha:\n");
      for(j = 0; j < Qmd; j++)
	printf("%lf ", logalpha(t, j));
      printf("\n");
    }
  }  // close the forward induction loop   
  mpfr_clear(var11);
  mpfr_clear(var21);
  mpfr_clear(summation3);
  
  /* forward termination */
  double ll = 0; // total log likelihood of all observation given this HMM
  for(i = 0; i < Qmd; i++){
    ll += logalpha(T-1, i);
  }
  ///===================================================================
  // for(i = 0; i < 100; i++){
  //   for(j = 0; j < Qmd; j++)
  //     printf("%lf ", logalpha(i, j));
  //   printf("\n");
  // }
  printf("\nprinting last column of logalpha...\n");
  for(i = 1; i < 6; i++){
    for(j = 0; j < Qmd; j++)
      printf("%lf ", logalpha(T-i, j));
    printf("\n");
  }
  printf("total loglikelihood: %lf\n", ll);
  ///===================================================================  
  double sum = 0;
  /* calculate logalpha last row sum */  
  for(i = 0; i < Qmd; i++)
    sum += logalpha(T-1, i);
  ll = sum;
  printf("LL: %lf........\n", ll);
  
  /* backward initilization */
  /// intialize mpfr variables  
  mpfr_t summation;
  mpfr_init(summation);
  mpfr_t var1, var2;
  mpfr_init(var1);
  mpfr_init(var2);
  mpfr_set_d(summation, 0.0, MPFR_RNDN);
  mpfr_set_d(var1, 0.0, MPFR_RNDN);
  mpfr_set_d(var2, 0.0, MPFR_RNDN);
  
  printf("backward initialization...\n");
  mpfr_set_d(summation, 0.0, MPFR_RNDN);
  double *temp = (double *)calloc(Qmd, sizeof(double ));
  
  for(i = 0; i < Qmd; i++)
    temp[i] = logalpha(T-1, i);
  
  for(i = 0; i < Qmd; i++){
    //double elem = logalpha(T-1, i);
    double elem = temp[i-1];
    mpfr_set_d(var2, elem, MPFR_RNDN);
    mpfr_exp(var1, var2, MPFR_RNDN);
    mpfr_add(summation, summation, var1, MPFR_RNDN);
  }
  // take logarithm
  mpfr_log(summation, summation, MPFR_RNDN);
  double sum2 = mpfr_get_d(summation, MPFR_RNDN);
  for(i = 0; i < Qmd; i++){    
    logg(T-1, i) = logalpha(T-1, i) - sum2 ;
  }

  // gamma matrix
  for(j = 0; j < Q; j++){
    gamma(j, T-1) = exp(logp_k(j, T-1) + logg(T-1, j));
  }
  mat lognewxi(Qmd, Qmd); // declare lognewxi matrix 
  /* backward induction */
  printf("backward induction in progress...\n");
  for(t = T-2; t >= 0 ; t--){
    for(j = 0; j < Qmd; j++){
      vec v1(Qmd);
      vec v2(Qmd);
      vec v3(Qmd);
      sum = 0;
      for(i = 0; i < Qmd; i++)
	v1(i) = logA(j, i);
      for(i = 0; i < Qmd; i++)
	v2(i) = logbeta(t+1, i);
      for(i = 0; i < Qmd; i++)
	v3(i) = new_logp(t+1, i);
      // add all three vectors
      for(i = 0; i < Qmd; i++)
	v1(i) += v2(i) + v3(i);
      mpfr_set_d(summation, 0.0, MPFR_RNDN);
      for(i = 0; i < Qmd; i++){
	double elem = v1(i);
	mpfr_set_d(var2, elem, MPFR_RNDN);
	mpfr_exp(var1, var2, MPFR_RNDN);
	mpfr_add(summation, summation, var1, MPFR_RNDN);
      }      
      mpfr_log(summation, summation, MPFR_RNDN);
      sum2 = mpfr_get_d(summation, MPFR_RNDN);
      logbeta(t, j) = sum2;
    }

    // computation of log(gamma) is now possible called logg here
    for(i = 0; i < Qmd; i++){
      logg(t, i) = logalpha(t, i) + logbeta(t, i);
    }
    mpfr_set_d(summation, 0.0, MPFR_RNDN);
    for(i = 0; i < Qmd; i++){
      double elem = logg(t, i);
      mpfr_set_d(var2, elem, MPFR_RNDN);
      mpfr_exp(var1, var2, MPFR_RNDN);
      mpfr_add(summation, summation, var1, MPFR_RNDN);
    }
    mpfr_log(summation, summation, MPFR_RNDN);
    sum2 = mpfr_get_d(summation, MPFR_RNDN);
    for(i = 0; i < Qmd; i++)
      logg(t, i) = logg(t, i) - sum2;
    
    // finally the gamma_k is computed (called gamma here )
    mpfr_set_d(summation, 0.0, MPFR_RNDN);
    for(j = 0; j < Q; j++){
      // for(i = j*MIN_DUR; i < (j+1) * MIN_DUR; i++){
      // 	sum += exp(logg(t, i));
      // }
      gamma(j, t) = exp( logp_k(j, t) + logg(t, j) );
    }    
    /* for the EM algorithm we need the sum
       over xi all over t */
    // replicate logalpha(t, :)' matrix along columns
    mat m1(Qmd, Qmd);
    for(i = 0; i < Qmd; i++){
      for(j = 0; j < Qmd; j++){
	m1(i, j) = logalpha(t, i);
      }
    }
    // replicate logbeta matrix
    vec v1(Qmd);
    for(i = 0; i < Qmd; i++)
      v1(i) = logbeta(t+1, i);
    vec v2(Qmd);
    for(i = 0; i < Qmd; i++)
      v2(i) = new_logp(t+1, i);
    vec v3(Qmd);
    for(i = 0; i < Qmd; i++)
      v3(i) = v1(i) + v2(i);
    // replicate v3 row vector along all rows of matrix m2
    mat m2(Qmd, Qmd);
    for(i = 0; i < Qmd; i++){
      for(j = 0; j < Qmd; j++){
	m2(i, j) = v3(i);
      }
    }
    
    // add both matrices m1 and m2 
    mat m3(Qmd, Qmd);
    m3 = m1 + m2;  // can do direct addition   
    ///mat lognewxi(Qmd, Qmd); // declare lognewxi matrix 
    lognewxi.zeros();
    lognewxi = m3 + logA;
    // add new sum to older sumxi
    /// first subtract total sum from lognewxi
    mpfr_set_d(summation, 0.0, MPFR_RNDN);
    for(i = 0; i < Qmd; i++){
      for(j = 0; j < Qmd; j++){
	double elem = lognewxi(i, j);
	mpfr_set_d(var2, elem, MPFR_RNDN);
	mpfr_exp(var1, var2, MPFR_RNDN);
	mpfr_add(summation, summation, var1, MPFR_RNDN);	
	//sum += exp(lognewxi(i, j));
      }
    }
    // now take the logarithm of sum
    mpfr_log(summation, summation, MPFR_RNDN);
    sum2 = mpfr_get_d(summation, MPFR_RNDN);
    // subtract sum from lognewxi
    for(i = 0; i < Qmd; i++){
      for(j = 0; j < Qmd; j++){
	lognewxi(i, j) = lognewxi(i, j) - sum2;
      }
    }
    mat newxi(Qmd, Qmd);
    newxi = lognewxi;
    // add sumxi and newlogsumxi
    /// take exponential of each element
    for(i = 0; i < Qmd; i++){
      for(j = 0; j < Qmd; j++){
	newxi(i, j) = exp(newxi(i, j));
      }
    }
    sumxi = sumxi + newxi;
  } // close the backward induction loop 
  
  /* handle annoying numerics */
  /// calculate sum of lognewxi along each row (lognewxi is already modified in our case)
  
  for(i = 0; i < Qmd; i++){
    mpfr_set_d(summation, 0.0, MPFR_RNDN);
    for(j = 0; j < Qmd; j++){
      //sum += lognewxi(i, j);
      double elem = lognewxi(i, j);
      mpfr_set_d(var2, elem, MPFR_RNDN);
      mpfr_exp(var1, var2, MPFR_RNDN);
      mpfr_add(summation, summation, var1, MPFR_RNDN);
    }
    sum2 = mpfr_get_d(summation, MPFR_RNDN);
    gamma1(i) = sum2;
  }
  // normalize gamma1 which is prior and normalize sumxi which is transition matrix
  sum = 0;
  for(i = 0; i < Qmd; i++)
    sum += gamma1(i);
  for(i = 0; i < Qmd; i++)
    gamma1(i) /= sum;  
  // transition probability matrix will be normalized in train_hmm function
  
  /// clear mpfr variables
  mpfr_clear(summation);
  mpfr_clear(var1);
  mpfr_clear(var2);
  printf("forward-backward algorithm calculation is done...\n");
  /* finished forward-backward algorithm */    
  return ll;
}
