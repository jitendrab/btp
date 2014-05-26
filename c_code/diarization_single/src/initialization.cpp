/*
Written By: Jitendra Prasad Keer
BTech, CSE, IIT Mandi
 */


/******************************************************************************
   InitializeGMMs : Uniform Initialization of all GMMs
   inputs : pointer to all feature vectors(VECTORS_OF_F_VECTORS), dimension of feature vectors, 
   totalNumFeatures, numStates
   outputs : Initialized GMMs with means vector and variance vector, posterior probabilities 
   all feature vectors
******************************************************************************/
int InitializeGMMs(VECTOR_OF_F_VECTORS *features, int Dim, int totalNumFeatures, int *numStates){
  
  double                                **posterior;
  int                                  i = 0, j = 0;
  int *path          = (int *)calloc(totalNumFeatures, sizeof(int ));
  //=====================================================================================
  //allocate memory  
  posterior = (double **)calloc(*numStates, sizeof(double *));
  for(i = 0; i < (*numStates); i++){
    posterior[i] = (double *)calloc(totalNumFeatures, sizeof(double));
  }
  for(i = 0; i < (*numStates); i++)
    for(j = 0; j < totalNumFeatures; j++)
      posterior[i][j] = 0.0;
  
  for(i = 0; i < totalNumFeatures; i++)
    features[i]->numElements = DIM;
  //=======================================================================================
  //  Initialize Minimum duration hmm
  
  // number of rows and columns in transition matrix
  int rowsPrior = (*numStates) * MIN_DUR;
  //number of columns and rows in prior matrix
  int rowsTrans = (*numStates) * MIN_DUR;
  int colsTrans = (*numStates) * MIN_DUR;

  struct timespec requestStart, requestEnd;
  clock_gettime(CLOCK_REALTIME, &requestStart);

  ESHMM *mdHMM                      = new ESHMM(*numStates, DIM, rowsTrans, colsTrans, rowsPrior); // calling constructor

  clock_gettime(CLOCK_REALTIME, &requestEnd);
  double accum = (requestEnd.tv_sec - requestStart.tv_sec) + (requestEnd.tv_nsec - requestStart.tv_nsec) / BILLION ;
  printf("time taken in constructor: %3.5lf\n", accum);
  
  printf("no of states: %d\n", mdHMM->hmmStates);
  mdHMM->hmmStates = *numStates;
  // Build GMM for each state  
  int numFeatures = totalNumFeatures/(*numStates); //last state should have all remaining features
  mat                    tempFeatures(numFeatures, DIM);
  tempFeatures.zeros();
  for(i = 0; i < (*numStates) ; i++){
    printf("\n performing uniform Initialization of GMM: %d....\n", i);
    int d = 0;
    mdHMM->numElemEachState[i] = numFeatures; // initial number of elements in each state
        
    for(j = 0; j < numFeatures; j++){
      for(d = 0; d < DIM; d++){
	tempFeatures(j,d) = features[i*numFeatures + j]->array[d];
      }
      path[i*numFeatures + j] = i;
    }
    clock_gettime(CLOCK_REALTIME, &requestEnd);
    
    // Estimate the Probability density funtion
    mdHMM->HMMstates[i].Estimate(tempFeatures.t());

    clock_gettime(CLOCK_REALTIME, &requestEnd);
    accum = (requestEnd.tv_sec - requestStart.tv_sec) + (requestEnd.tv_nsec - requestStart.tv_nsec) / BILLION ;
    printf("time in Estimation, %d state: %3.5lf\n", i, accum);

    //Cov = mdHMM->HMMstates[i].Covariance();    
    //Cov.print("printing Cov:\n");   
  }  
  printf("Initialization complete\n\n");  
  // Basic HMM is Complete
  //printHMM(mdHMM);
  // convert to minimum duration hmm
  hmm2MinDurationHMM(mdHMM);
  //mdHMM->printPriorMat(mdHMM);
  mdHMM->printTransMat(mdHMM);
  /// extract modifiable copy of Mean
  arma::vec &Mean = mdHMM->HMMstates[0].Mean();
  Mean.print("Printing Mean:\n");
  //for(i = 0; i < DIM; i++)
  // m(i) = 5.0;
  /// print new Mean
  
  vec x(DIM);
  x.zeros();
  for(j = 0; j < DIM; j++)
    x(j) = features[0]->array[j];
  double prob = mdHMM->HMMstates[2].Probability(x) ;
  printf("prob: %lf\n", prob);
  prob = log(prob);
  printf("p: %lf\n", prob);

  double neg_inf = -INFINITY;
  //  assert(std::isinf(neg_inf));
  
  ComputePosteriorProb(mdHMM, features, totalNumFeatures, posterior);
  mdHMM->hmmLogViterbiWithMinDur(path, mdHMM, totalNumFeatures, posterior);
  
  printf("neg_inf: %lf\n", neg_inf);
  ClusteringAndMerging(features, path, mdHMM, totalNumFeatures, posterior);
  // print the path obtained  
  printf("Everything is done......\n");
  return 0;
}
/*******************************************************************************/
void ESHMM::FindNumElemEachState(ESHMM *mdHMM, int *path, int totalNumFeatures){
  int i = 0;
  int numStates = mdHMM->hmmStates;
  for(i = 0; i < numStates; i++)
    mdHMM->numElemEachState[i] = 0;
  
  for(i = 0; i < totalNumFeatures; i++){
    int idx = path[i];
    mdHMM->numElemEachState[idx]++;
  }
  for(i = 0; i < numStates; i++){
    printf("features in state, %d : %d\n", i, mdHMM->numElemEachState[i]);
  }
  return;
}
/*******************************************************************************/
void ClusteringAndMerging(VECTOR_OF_F_VECTORS *features, int *path, ESHMM *mdHMM, int totalNumFeatures, double **posterior){
  int flag = 1, i = 0, j = 0;
  int s = 0;
  double **deltaBIC = (double **)calloc(MAX_NUM_STATES, sizeof(double *));
  for(i = 0; i < MAX_NUM_STATES; i++){
    deltaBIC[i] = (double *)calloc(MAX_NUM_STATES, sizeof(double ));
  }
  
  while(flag && mdHMM->hmmStates > 1){
    // Merge two States, drop any empty state, calculate BIC and find best match
    /* decrease number of states */
    int dropped_state = hmmMergeTwoStates(mdHMM, path, totalNumFeatures, deltaBIC, features);
    
    if(dropped_state == -1 || mdHMM->hmmStates < 2){
      printf("we have received your signal...stopping...\n");
      flag = 0;
      printf("final no of states: %d\n", mdHMM->hmmStates);
    }else{
      //create a new HMM and train it using Baum Welch algorithm, with minimum duration
      printf("creating a new MDHMM ...\n");
      // creating a new HMM
      // drop states which has less number of elements
      int empty = 1;
      while(empty){
	mdHMM->FindNumElemEachState(mdHMM, path, totalNumFeatures);
	for(s = 0; s < mdHMM->hmmStates; s++){
	  if(mdHMM->numElemEachState[s] == 0){
	    mdHMM->DropEmptyStates(mdHMM, path, totalNumFeatures);	    
	  }
	}
	empty = 0;
	for(s = 0; s < mdHMM->hmmStates; s++){
	  if(mdHMM->numElemEachState[s] == 0)
	    empty = 1;
	}
      }
      // modify current HMM, remove one state, update prior and transition matrices
      printf("modifying current HMM, decreasing no of states...\n");
      // pop one state from the state vector
      mdHMM->HMMstates.pop_back();
      /* convert modified HMM to Min-Duration HMM */
      hmm2MinDurationHMM(mdHMM);
      // Initialize HMM using path 
      mdHMM->InitializeHMMusingPath(mdHMM, path, features, totalNumFeatures);
      ComputePosteriorProb(mdHMM, features, totalNumFeatures, posterior); // log probs
      // train new Min duration HMM using baum welch algorithm
      mdHMM->trainMDHMM(mdHMM, features, totalNumFeatures, posterior);
      // calculate new optimal path using viterbi algorithm
      mdHMM->hmmLogViterbiWithMinDur(path, mdHMM, totalNumFeatures, posterior);
      empty = 1;
      while(empty){
	mdHMM->FindNumElemEachState(mdHMM, path, totalNumFeatures);
	for(s = 0; s < mdHMM->hmmStates; s++){
	  if(mdHMM->numElemEachState[s] == 0){
	    mdHMM->DropEmptyStates(mdHMM, path, totalNumFeatures);
	  }
	}
	empty = 0;
	for(s = 0; s < mdHMM->hmmStates; s++){
	  if(mdHMM->numElemEachState[s] == 0)
	    empty = 1;
	}
      }
      
    } // close else
  } // close while 
  for(i = 0; i < totalNumFeatures; i++)
    printf("%d ", path[i]);  
  printf("\n\n");
  // write RTTM file
  return;
}
/*******************************************************************************
drop any state which contain no of elements less than MIN_DUR

*******************************************************************************/
void ESHMM::DropEmptyStates(ESHMM *mdHMM, int *path, int T){
  mdHMM->FindNumElemEachState(mdHMM, path, T);
  int states = mdHMM->hmmStates;
  int s = 0, i = 0;
  for(s = 0; s < states; s++){
    // if number of elements are less than MIN_DUR
    if(mdHMM->numElemEachState[s] < MIN_DUR){
      // consider dropping this state
      printf("dropping Empty state: %d....\n", s);
      for(i = 0; i < T; i++){
	if(path[i] > s)
	  path[i] = path[i] - 1;
      }
      // decrease number of states
      mdHMM->hmmStates -= 1;
      break;
    }
  }
  // modify prior 
  mdHMM->rowsPrior = mdHMM->hmmStates * MIN_DUR;
  // modify transition matrix
  mdHMM->rowsTrans = mdHMM->colsTrans = mdHMM->hmmStates * MIN_DUR;
}
/*******************************************************************************/
void ESHMM::InitializeHMMusingPath(ESHMM *mdHMM, int *path, VECTOR_OF_F_VECTORS *features, int totalNumFeatures){
  printf("initializing new MDHMM using new viterbi state path...\n");
  int states = mdHMM->hmmStates;
  mdHMM->FindNumElemEachState(mdHMM, path, totalNumFeatures);
  int i = 0, j = 0, s = 0, k = 0;
  for(s = 0; s < states; s++){
    int numFeatures = mdHMM->numElemEachState[s];
    mat tempFeatures(numFeatures, DIM);
    int count = 0;
    for(k = 0; k < totalNumFeatures; k++){
      if(path[k] == s){
	int d = 0;
	for(d = 0; d < DIM; d++)
	  tempFeatures(count, d) = features[k]->array[d];
	count++;
      }
    }
    // Estimate the Probability density funtion
    mdHMM->HMMstates[s].Estimate(tempFeatures.t());
  }
  //change prior 
  mdHMM->rowsPrior = states * MIN_DUR;
  // change transition matrix
  mdHMM->rowsTrans = states * MIN_DUR;
  mdHMM->colsTrans = states * MIN_DUR;
}
/*******************************************************************************/
void ESHMM::trainMDHMM(ESHMM *mdHMM, VECTOR_OF_F_VECTORS *features, int totalNumFeatures, double **posterior){
  printf("training min-duration hmm using all features...\n");
  int i = 0, j = 0;
  // you can add stopping criterion as fun argument which includes maxIters and minLLimprovement
  // in each iteration
  // Default Values: 
  double minLLimpr = 1e-5;
  int maxIter = 5;
  double reg = 0.001;
  printf("default parameters used:\n minLLimpr: %lf \n maxIterations: %d\n reg: %lf\n", minLLimpr, maxIter, reg);
  // By default we are dealing with minimum duration HMMs
  int Q = mdHMM->hmmStates; 
  int Qmd = Q * MIN_DUR;
  /* initialization */
  int T = totalNumFeatures; // DIM is already defined  
  int iter = 0;
  double LL2 = -INFINITY;
  /* matrices to be passed as arguments to forward-backward */  
  mat sumxi(Qmd, Qmd); // make life simpler later on 
  mat gamma(Q, T);  // 
  /* create new gamma1 matrix */    
  rowvec gamma1(Qmd);
  
  while(iter < maxIter){
    printf("hmm iter:%d \n", iter);
    /* compute log sum probabilities using Min-Duration constraint*/  /* */
    ComputePosteriorProb(mdHMM, features, totalNumFeatures, posterior);
    /* forward backward */
    double LL1 = LL2;
    sumxi.zeros();
    gamma.zeros();
    gamma1.zeros();
    LL2 = mdHMMLogForwardBackward(mdHMM, features, posterior, totalNumFeatures, gamma, gamma1, sumxi);
    printf("HMM iteration: %d   LL2: %lf  LL1: %lf\n", iter, LL2, LL1);    
    /* is it going into the right direction */ 
    if(LL2 < LL1){
      printf("LL difference: %lf and minLLimpr: %lf \n", LL2 - LL1, minLLimpr);
      printf("LLimpr is less than minLLimpr...stopping...\n");
      break;
    }
    /* estimate the model parameters again: */
    /* - prior and transition probabilities */    
    /// mdHMM->prior = gamma1'; // gamma1 is a row vector 
    for(i = 0; i < mdHMM->rowsPrior; i++)
      mdHMM->prior->array[i] = log(gamma1(i));
    /// mdHMM->trans = sumxi.t();
    /* normalise transition probability matrix */
    sumxi = sumxi.t(); // take the transpose    
    double sum = 0.0;
    for(i = 0; i < Qmd; i++){
      sum = 0;
      for(j = 0; j < Qmd; j++)
	sum += sumxi(i, j);
      //now divide each element in row with sum
      for(j = 0; j < Qmd; j++)
	sumxi(i, j) = sumxi(i, j) / sum;
    }
    gamma1.print("printing gamma1:\n");
    sumxi.print("printing sumxi:\n");
    // assign sumxi to transition matrix
    for(i = 0; i < Qmd; i++){
      for(j = 0; j < Qmd; j++){
	mdHMM->trans[i]->array[j] = log(sumxi(i, j));
      }
    }
    
    /*  state models, (prior, mean, covariance matrix)*/ 
    mat X(T, DIM);
    
    for(i = 0; i < Q; i++){
      int mi = 1;
      double Z = 0;
      /// sum the gamma along column
      for(j = 0; j < T; j++)
	Z += gamma(i, j);
      /// now make a matrix biggamma
      /// replicate gamma into each dimension
      if(Z == 0) /// if zero
	Z == 1;
      
      ///model.pdf{i}.prior  = 1
      int d = 0;
      for(j = 0; j < T; j++){
	for(d = 0; d < DIM; d++){
	  X(j, d) = ( features[j]->array[d] * gamma(i, j) ) / Z;
	}
      }
      /// now to obtain new mean, sum along dim1
      vec temp_mean(DIM);
      for(d = 0; d < DIM; d++){
	sum = 0.0;
	for(j = 0; j < T; j++){
	  sum += X(j, d);
	}
	temp_mean(d) = sum;
      }
      temp_mean.print("printing new mean: \n");
      // extract new modfiable copy of Mean and modify the mean
      arma::vec &Mean = mdHMM->HMMstates[i].Mean();
      /// modify the Mean
      for(j = 0; j < DIM; j++)
	Mean(j) = temp_mean(j);
      printf("New mean has been assigned to the state: %d\n", i);
      /* obtain new covariance matrix */
      mat dx(T, DIM);
      for(j = 0; j < T; j++){
	for(d = 0; d < DIM; d++){
	  dx(j, d) = features[j]->array[d] - temp_mean(d);
	}
      }
      mat dx2 = dx;
      /// multiply replicated gamma and dx2
      for(j = 0; j < T; j++){
	for(d = 0; d < DIM; d++)
	  dx2(j, d) = dx2(j, d) * gamma(i, j);
      }
      mat newcov = (dx.t() * dx2 ) / Z + reg * eye(DIM, DIM) ;
      /// extract modifiable covariance matrix from vector and assign new covariance matrix
      arma::mat &Cov = mdHMM->HMMstates[i].Covariance();
      /// modify new covariance matrix
      for(j = 0; j < DIM; j++){
	for(d = 0; d < DIM; d++)
	  Cov(j, d) = newcov(j, d);
      }
      printf("New Covariace has been assigned to the state: %d\n", i);
      /// its done
    }
    /// increment no of iterations
    iter++;
  }
}


/******************************************************************************
Calculate deltaBIC between each pair of states, find best pair of states to merge, 
Merge two states, change state in path, decrease number of states by one 
******************************************************************************/
int hmmMergeTwoStates(ESHMM *mdHMM, int *path, int totalNumFeatures, double **deltaBIC, VECTOR_OF_F_VECTORS *features){
  int i = 0, j = 0;
  int states = mdHMM->hmmStates;  
  printf("calculating delta BIC to select best pair to merge...\n");
  printf("no of states: %d\n", states);
  // drop any state which is empty, not assigned a single feature or element
  mdHMM->FindNumElemEachState(mdHMM, path, totalNumFeatures);
  for(i = 0; i < states; i++){
    if(mdHMM->numElemEachState[i] <= MIN_DUR)
      printf("no of elems in %d: %d... consider dropping the state...\n", i, mdHMM->numElemEachState[i]);
  }

  // calculate delta BIC for each pair of state
  CalculateDeltaBIC(mdHMM, deltaBIC, features, totalNumFeatures, path);
  // find the minimum delta bic pair of state
  int min_i = 0, min_j = 1;
  double min = deltaBIC[0][1];
  for(i = 0; i < states; i++){
    for(j = i+1; j < states; j++){
      if(deltaBIC[i][j] <= min){
	min = deltaBIC[i][j];
	min_i = i;
	min_j = j;
      }
    }
  }
  printf("min deltaBIC: %lf, min_i: %d  min_j: %d\n", min, min_i, min_j);
  if(deltaBIC[min_i][min_j] > 0){
    printf("we need to stop now...merging is not possible...\n");
    return -1;
  }else{
    // merge state_i and state_j
    int numFeatures = mdHMM->numElemEachState[min_i] + mdHMM->numElemEachState[min_j];
    mat tempFeatures(numFeatures, DIM);
    for(i = 0; i < totalNumFeatures; i++){
      if(path[i] == min_j){
	path[i] = min_i;
      }
      else if(path[i] > min_i){
	path[i] = path[i] - 1;
      }
    }
    // change number of elements in each state
    mdHMM->numElemEachState[min_i] += mdHMM->numElemEachState[min_j];
    for(i = min_j + 1; i < mdHMM->hmmStates - 1; i++){
      mdHMM->numElemEachState[i] = mdHMM->numElemEachState[i+1];
    }
    // Build a New HMM according to new path, Initialize and train the new HMM
    // decrease HMM states by one
    mdHMM->hmmStates -= 1;
    return 1;
  }    
}

/*******************************************************************************/
void CalculateDeltaBIC(ESHMM *mdHMM, double **deltaBIC, VECTOR_OF_F_VECTORS *features, int T, int *path){
  int i = 0, j = 0, s = 0, k = 0;
  int hmmStates = mdHMM->hmmStates;
  double *det = (double *)calloc(hmmStates, sizeof(double ));
  double L0 = 0, L1 = 0;
  for(s = 0; s < hmmStates; s++){
    //calculate determinant of covariance matrix for each state
    mat &Cov = mdHMM->HMMstates[s].Covariance();
    double sign;
    double value = 0;
    log_det(value, sign, Cov);
    det[s] = value;
  }

  for(i = 0; i < hmmStates; i++){
    for(j = i+1; j < hmmStates; j++){
      //model state_i and state_j using single gaussian
      int numFeatures = mdHMM->numElemEachState[i] + mdHMM->numElemEachState[j];
      mat tempFeatures(numFeatures, DIM);
      int count = 0;
      for(k = 0; k < T; k++){
	if(path[k] == i || path[k] == j){
	  int d = 0;
	  for(d = 0; d < DIM; d++)
	    tempFeatures(count, d) = features[k]->array[d];
	  count++;
	}
      }
      
      GaussianDistribution gaussian (DIM);
      gaussian.Estimate(tempFeatures.t());
      // find covariance matrix
      mat &Cov = gaussian.Covariance();
      //Cov.print();
      // calculate det of covariance matrix
      double detSigma = 0;
      double sign;
      log_det(detSigma, sign, Cov);
      int n_i = mdHMM->numElemEachState[i];
      int n_j = mdHMM->numElemEachState[j];
      double penalty = LAMBDA * 0.5 * (DIM + 0.5 * DIM * (DIM + 1)) * log(n_i + n_j);
      printf("det[%d]: %lf det[%d]: %lf detSigma: %lf \n", i, det[i], j, det[j], detSigma);
      deltaBIC[i][j] = (n_i + n_j) * detSigma - n_i * det[i] - n_j * det[j] - penalty;
      printf("i:%d  j:%d  deltaBIC: %lf P: %lf, n_i:%d  n_j:%d\n\n", i, j , deltaBIC[i][j], penalty, n_i, n_j);
    }
  }
  printf("delta bic calculation is completed...\n");
}
/*******************************************************************************/
void ESHMM::hmmLogViterbiWithMinDur(int *path, ESHMM *mdHMM, int T, double **B){
  int              Q = mdHMM->hmmStates * MIN_DUR;
  mat              delta(T, Q);
  mat              psi(T, Q);
  mat              scale(T, 1);
  int              i = 0, j = 0, s = 0, t = 0, k = 0;
  int              numStates = mdHMM->hmmStates;
  // minimum duration viterbi hence modify B(posterior) prob matrix
  mat              new_B(T, Q*MIN_DUR);
  
  for(i = 0; i < numStates; i++){
    for(j = 0; j < T; j++){
      for(k = i*MIN_DUR; k < (i+1)*MIN_DUR; k++){
	new_B(j, k) = B[i][j];
      }
    }
  }
  
  // first step
  /// prior and transition probabilities are in log so we will use log from now onwards
  /* Using logarithm for exact numerical calculations  */
  for(i = 0; i < Q; i++){
    delta(0, i) = mdHMM->prior->array[i] + new_B(0, i);
  }
  printf("calculating delta and psi matrix:...\n");
  for(t = 1; t < T; t++){
    for(j = 0; j < Q; j++){
      double max = -INFINITY; // just initializing max 
      int max_s = 0;
      for(i = 0; i < Q; i++){
	if((delta(t-1, i) + mdHMM->trans[i]->array[j]) > max){
	  max = delta(t-1, i) + mdHMM->trans[i]->array[j];
	  max_s = i;
	}
      }
      delta(t, j) = max + new_B(t, j);
      psi(t, j) = max_s;
    }
    
    // you can check airthmetic integrity of calculation becasue of very low values
  }
  printf("finding the path backward...\n");
  double max = delta(T-1, 0);
  int max_s = 0;
  for(j = 1; j < Q; j++){
    if(delta(T-1, j) > max){
      max = delta(T-1, j);
      max_s = j;
    }
  }
  path[T-1] = max_s;
  printf("path[T-1]: %d\n", path[T-1]);
  for(t = T-2; t >= 0; t--){
    int idx = 0;
    idx = path[t+1];
    //printf("idx:%d \n", idx);
    path[t] = psi(t+1, idx);
  }
  
  // condense the path for minimum duration
  for(i = 0; i < T; i++){
    path[i] = ceil(path[i] / MIN_DUR);
  }
  printf("viterbi path finding is complete...\n");
  return;
}
/******************************************************************************* 
   ComputePosteriorProb() : compute posterior probabilities (log probabilities)
   features - pointer to posterior matrix , pointer to all feature vectors, MDHMM, totalNumFeatures, 
   posterior - posterior probability matrix
   outputs : compute posterior probability for each feature vector 
******************************************************************************/
void ComputePosteriorProb( ESHMM *mdHMM, VECTOR_OF_F_VECTORS *features, int totalNumFeatures, double **posterior){
  vec               x(DIM);
  int               k  = 0, mixCount = 0, i = 0, j = 0, s = 0;
  int hmmStates     = mdHMM->hmmStates;
  double            p = 0.0;
  //compute posterior probabilities
  printf("calculating posterior probabilities for all feature vectors.......\n");
  printf("no of hmm states: %d\n", hmmStates);
  for(i = 0; i < totalNumFeatures; i++){
    for(j = 0; j < hmmStates; j++){
      posterior[j][i] = 0.0;
    }
  }
  /// testing for singularity, 
  /// print det cov for each state
  for(s = 0; s < hmmStates; s++){
    mat Cov = mdHMM->HMMstates[s].Covariance();
    double sign, det = 0;
    log_det(det, sign, Cov);
    printf("log_det(%d): %lf\n", s, det);
  }
  
  for(s = 0;s < (hmmStates); s++){
    printf("\n s: %d\n", s);
    for(i = 0; i < totalNumFeatures; i++){
      //printf("%d ", i);
      for(k = 0; k < DIM; k++)
	x(k) = features[i]->array[k];
      //x.print("x is:\n");
      //printf("\nfor feature vector: %d\n", i);
      p = mdHMM->HMMstates[s].Probability(x);
      p = log(p);
      posterior[s][i] = p;
    }
    //printf("%f  ", posterior[j][i]);
  }
  //printf("\n");
  printf("finished...calculation of posterior probabilities.\n");
}
/******************************************************************************* */
void ESHMM::printPriorMat(ESHMM *mdHMM){
  int i = 0;
  int rows = mdHMM->prior->numElements;
  printf("Printing Prior:\n");
  for(i = 0; i < rows; i++)
    printf("%lf ", mdHMM->prior->array[i]);
  printf("\n");
}
/******************************************************************************* */
void ESHMM::printTransMat(ESHMM *mdHMM){
  printf("printing Transition Matrix:\n");
  int rows = mdHMM->rowsTrans;
  int cols = mdHMM->colsTrans;
  int i = 0, j= 0;
  for(i = 0; i < rows; i++){
    for(j = 0; j < cols; j++){
      printf("%lf   ", mdHMM->trans[i]->array[j]);
    }
    printf("\n");
  }
  printf("\n");
}
/******************************************************************************* */
void hmm2MinDurationHMM(ESHMM *mdHMM){
  printf("converting HMM to min duration hmm\n");
  int i = 0, j = 0;      
  int hmmStates = mdHMM->hmmStates; // very important to correctly identify no of states
  vec prior(hmmStates);
  double prob = (double)1.0/hmmStates;  /// we will use log probability 
  prior.fill(prob);

  mat trans(hmmStates, hmmStates);
  trans.fill(prob);
  trans.print("trans:\n");
  // make all zeros in transition matrix
  for(i = 0; i < hmmStates * MIN_DUR; i++){
    for(j = 0; j < hmmStates * MIN_DUR; j++)
      mdHMM->trans[i]->array[j] = 0.0;
  }
  /// make all zeros in prior matrix
  for(i = 0; i < hmmStates * MIN_DUR; i++)
    mdHMM->prior->array[i] = 0.0;
  
  // create prior matrix  
  printf("states: %d  prior: %f\n", hmmStates, prob);
  for(i = 0; i < hmmStates*MIN_DUR; i += MIN_DUR){
    mdHMM->prior->array[i] = prob;
  }
  /// take logarithm of prior probabilities
  for(i = 0; i < hmmStates * MIN_DUR; i++){
    if(mdHMM->prior->array[i] == 0.0)
      mdHMM->prior->array[i] = -INFINITY;
    else
      mdHMM->prior->array[i] = log(mdHMM->prior->array[i]);
  }
  
  // create a topeliz matrix
  for(i = 0, j = 1; i < hmmStates*MIN_DUR, j < hmmStates*MIN_DUR; i++, j++){
    mdHMM->trans[i]->array[j] = 1;
  }
  
  // now copy the elements on the right spot
  for(i = 1; i <= hmmStates; i++){
    mdHMM->trans[i*MIN_DUR - 1]->array[i*MIN_DUR -1] = trans(i-1, i-1); //hmm.trans[i-1][i-1] all are equal to prob, might not be equal
    for(j = i+1; j <= hmmStates; j++){
      mdHMM->trans[(i)*MIN_DUR -1]->array[(j-1) * MIN_DUR ] = trans(i-1, j-1); // Warning change transition probabilties to real probabilities 
      mdHMM->trans[(j)*MIN_DUR -1]->array[(i-1)*MIN_DUR] = trans(j-1, i-1);
    }
  }
  // model.trans = sparse(trans)   
  /* now take the log of transition matrix, define log(0.0) as -infinity */
  for(i = 0; i < mdHMM->rowsTrans; i++){
    for(j = 0; j < mdHMM->colsTrans; j++){
      if(mdHMM->trans[i]->array[j] == 0.0){
	mdHMM->trans[i]->array[j] = -INFINITY;
      }
      else
	mdHMM->trans[i]->array[j] = log(mdHMM->trans[i]->array[j]);
    }
  }
    
  mdHMM->MD = 250;
  return;
}   
/******************************************************************************* */

