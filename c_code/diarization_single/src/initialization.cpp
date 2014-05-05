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
  
  float                                *posterior[*numStates];
  int                                  i = 0, j = 0;
  int                                  *numElemEachState;
  
  //=====================================================================================
  //allocate memory  
  
  for(i = 0; i < (*numStates); i++){
    posterior[i] = (float *)calloc(totalNumFeatures, sizeof(float));
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
    printf("\nperforming uniform Initialization of GMM: %d....\n", i);        
    int d = 0;
    mdHMM->numElemEachState[i] = numFeatures; // initial number of elements in each state
        
    for(j = 0; j < numFeatures; j++){
      for(d = 0; d < DIM; d++){
	tempFeatures(j,d) = features[i*numFeatures + j]->array[d];
      }
    }
    clock_gettime(CLOCK_REALTIME, &requestEnd);
    
    // Estimate the Probability density funtion
    mdHMM->HMMstates[i].Estimate(tempFeatures.t());

    clock_gettime(CLOCK_REALTIME, &requestEnd);
    accum = (requestEnd.tv_sec - requestStart.tv_sec) + (requestEnd.tv_nsec - requestStart.tv_nsec) / BILLION ;
    printf("time in Estimation, %d state: %3.5lf\n", i, accum);
    //Mean = mdHMM->HMMstates[i].Mean();
    //Cov = mdHMM->HMMstates[i].Covariance();    
    //Mean.print("Printing Mean:\n");
    //Cov.print("printing Cov:\n");
    //print temp features
    /* for(i = 0; i < numFeatures; i++){
       for(j = 0; j < Dim; j++)
       printf("%f  ", tempFeatures[i]->array[j]);
       printf("\n");
       }
    */  
  }  
  printf("Initialization complete\n\n");
  // Basic HMM is Complete
  printHMM(mdHMM);
  // convert to minimum duration hmm
  hmm2MinDurationHMM(mdHMM);
  //ComputePosteriorProb(features,    posterior,      allMixtureMeans,       allMixtureVars,     numStates,    totalNumFeatures);  
  
  
  //ClusteringAndMerging( features,   mdHMM,    totalNumFeatures,   posterior, numElemEachState,    Pi);
  
  /*int elems_in_states[MAX_NUM_STATES] = {0};
  //classify using posterior probabilities
  for(i = 0; i < totalNumFeatures; i++){
    float max = -999999.0;
    int max_s = 0;
    for(j = 0; j < *numStates; j++){
      if(posterior[j][i] > max){
	max = posterior[j][i];
	max_s = j;
      }
    }
    elems_in_states[max_s] += 1;
  }
  for(i = 0; i < *numStates; i++)
    printf("elem in state: %d     %d \n", i, elems_in_states[i]);
  */

  printf("Everything is done......\n"); 
  return 0;
}
/******************************************************************************* */
void hmm2MinDurationHMM(ESHMM *mdHMM){
  int i = 0, j = 0;      
  // create prior matrix
  int numStates = mdHMM->numStates;
  float prob = (double)1.0/numStates);
  printf("states: %d  prior: %f\n", numStates, prob);
  for(i = 0; i < numStates*MIN_DUR; i += MIN_DUR){
    mdHMM->prior->array[i] = prob;
  }
  
  // create a topeliz matrix
  for(i = 0, j = 1; i < numStates*MIN_DUR, j < numStates*MIN_DUR; i++, j++){
    trans[i]->array[j] = 1;
  }
  
  // now copy the elements on the right spot
  for(i = 1; i <= numStates; i++){
    trans[i*MIN_DUR-1]->array[j*MIN_DUR-1] = prob; //hmm.trans[i-1][i-1] all are equal to prob
    for(j = i+1; j <= *numStates; j++){
      trans[i*MIN_DUR -1]->array[(j-1) * MIN_DUR ] = prob; // all transition prob are same in our case
      trans[j*MIN_DUR -1]->array[(i-1)*MIN_DUR ] = prob;
    }
  }
  // model.trans = sparse(trans) 
  mdHMM->prior = prior;
  mdHMM->trans = trans;
  mdHMM->MD = 250;
  return mdHMM;
} 
  
