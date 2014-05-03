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

  VECTOR_OF_F_VECTORS                  *allMixtureMeans, *allMixtureVars;
  
  float                                *posterior[*numStates];
  int                                  i = 0, j = 0;
  float                                *mixtureElemCount[*numStates];
  int                                  VQIter = 1, GMMIter = 2;//Number of iterations  
  int                                  *numElemEachState;
  //=====================================================================================
  //allocate memory  
  featuresForClustering = new F_VECTOR *[totalNumFeatures];
  for(i = 0; i < totalNumFeatures; i++){
    //featuresForClustering[i] = (F_VECTOR *)AllocFVector(DIM);
    featuresForClustering[i] = new F_VECTOR(DIM);
  }
  numElemEachState           = (int *)calloc(*numStates, sizeof(int));
        
  mixtureMeans               = (VECTOR_OF_F_VECTORS *)calloc(MAX_NUM_MIX, sizeof(VECTOR_OF_F_VECTORS ));
  mixtureVars                = (VECTOR_OF_F_VECTORS *)calloc(MAX_NUM_MIX, sizeof(VECTOR_OF_F_VECTORS ));
  
  for(i = 0; i < INITIAL_NUM_MIX * (*numStates); i++){
    //mixtureMeans[i]                  = (F_VECTOR *)AllocFVector(Dim);
    mixtureMeans[i]                  = new F_VECTOR(DIM);
    mixtureMeans[i]->numElements     = DIM;
    //mixtureVars[i]                   = (F_VECTOR *)AllocFVector(Dim);
    mixtureVars[i]                   = new F_VECTOR(DIM);
    mixtureVars[i]->numElements      = DIM;
  }
  for(i = 0; i < (*numStates); i++){
    mixtureElemCount[i] = (float *)calloc(INITIAL_NUM_MIX, sizeof(float));
  }
  
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
  int colsPrior = 1;  
  //number of columns and rows in prior matrix
  int rowsTrans = (*numStates) * MIN_DUR;
  int colsTrans = (*numStates) * MIN_DUR;
  
  // ESHMM                      *mdHMM;
  ESHMM *mdHMM                      = new ESHMM(*numStates, DIM, rowsTrans, colsTrans, rowsPrior, colsPrior); // calling constructor
  // Build GMM for each state  
  for(i = 0; i < (*numStates); i++){    
    printf("\nperforming uniform Initialization of GMM: %d....\n", i);    
    int numFeatures = totalNumFeatures/(*numStates); //last state should have all remaining features
    int d = 0;
    numElemEachState[i] = numFeatures; // initial number of elements in each state
    for(j = 0; j < numFeatures; j++){
      for(d = 0; d < DIM; d++){
	featuresForClustering[j]->array[d] = features[i*numFeatures + j]->array[d];
	featuresForClustering[j]->numElements = DIM;
      }
    }
    //print temp features
    /* for(i = 0; i < numFeatures; i++){
       for(j = 0; j < Dim; j++)
       printf("%f  ", tempFeatures[i]->array[j]);
       printf("\n");
       }
    */
    //    ComputeGMM(featuresForClustering,            numFeatures,              mixtureMeans,             mixtureVars,                       mixtureElemCount[i],                INITIAL_NUM_MIX,          VQIter,           GMMIter,                       probScaleFactor,     ditherMean,                         varianceNormalize,        time(NULL));
    int mixCount = 0, k=0;
    //store current mean and variance into all means and variance        
    //mdHMM->HMMstates[i]->mean->array[k]
    for(k = 0; k < Dim; k++){
      mdHMM->HMMstates[i]->mean->array[k] = mixtureMeans[0]->array[k];
      mdHMM->HMMstates[i]->cov[k]->array[k] = mixtureVars[0]->array[k];
    }
    printf("Initialization complete\n\n");
  }
  // Basic HMM is Complete
  printHMM(mdHMM);
  // convert to minimum duration hmm
  //  mdHMM = hmm2MinDurationHMM(mdHMM, DIM, numStates, MIN_DUR);
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

  printf("Everything is done...now SIGsegv will occur....\n");
  // FREE THE MEMORY
  //free(numElemEachState);
  //free(Pi);
  //free(numMixEachState);
  return 0;
}
/******************************************************************************* */
