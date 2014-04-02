#include "../include/preprocessing.h"

int main(int argc, char *argv[]){
  if(argc < 3){
    fprintf(stderr, "Usage: <binary> <path to feature file(HTK format)>  <path to .scp file>\n");
    exit(0);
  }
  // Variable declaration
  char                                 buffer[512];
  FILE                                 *feat_file, *scp_file;
  int                                  i = 0, j= 0, count = 0, TotalFeatures = 0,ret=0;
  char                                 line[512]; //to store the whole line    
  short int                            flag[MAX_NUM_FEATURES] = {0};
  VECTOR_OF_F_VECTORS                  *features; //check memory access (need to allocate more space)
  F_VECTOR                             *feat;
  int                                  States = 6;
  int                                  *numStates = &States;
  features = (VECTOR_OF_F_VECTORS *)calloc(MAX_NUM_FEATURES, sizeof(VECTOR_OF_F_VECTORS));
  feat = (F_VECTOR *)AllocFVector(DIM);
  for(i = 0; i < MAX_NUM_FEATURES; i++){
    features[i] = (F_VECTOR *)AllocFVector(DIM);
  }
  //open feature file(HTK format) and convert to ASCII and write it  
  snprintf(buffer, sizeof(buffer), "/usr/local/bin/HList -r %s > feature_file_all.txt", argv[1]);
  if(system(buffer) < 0){
    fprintf(stderr, "unable to open HTK binary file: HList or paths to files are wrong\nExiting\n");
    exit(0);
  }
  // open feature file
  feat_file = fopen("feature_file_all.txt", "r");
  if(feat_file == NULL){
    fprintf(stderr, "unable to open feature_file_all.txt\n");
    exit(0);
  }
  // open scp file
  scp_file = fopen(argv[2], "r");
  if(scp_file == NULL){
    fprintf(stderr, "unable to open scp file\n");
    exit(0);
  }    
  
  while((ret = fscanf(scp_file, "%s", line)) == 1){    
    if(ret == EOF)
      break;  
    char start[20], end[20];
    i = 0;
    while(line[i] != '_')
      i++;
    for(j = 0,i++; line[i] != '_'; j++,i++)
      start[j] = line[i];
    for(j=0, i++; line[i] != '='; i++, j++)
      end[j] = line[i];
    int s = atoi(start);
    int e = atoi(end);
    for(i = s; i <= e; i++)
      flag[i] = 1;    
    // printf("%d  %d \n", s, e);
  }
  
  //process feature file and drop all unvoiced frames         
  if(!features || !feat){
    fprintf(stderr, "unable to allocate memory\n");
    exit(0);
  }
  
  count = 0; //count actual feature vectors
  i = 0; // count file feature vectors
  while(!feof(feat_file) && !ferror(feat_file)){
    i++;    
    float number = 0;
    for(j = 0; j < DIM; j++){
      fscanf(feat_file, "%f", &number);
      feat->array[j] = number;
    }
    //check if i is voiced feature vector or not
    if(flag[i] == 0){
      continue;
    }
    else{
      //this feature vector is voiced 
      for(j = 0; j < DIM; j++)
	features[count]->array[j] = feat->array[j];
      features[count]->numElements = DIM;
      count++;
      //      printf("%d  ", i);
    }
  }
  TotalFeatures = count;
  printf("TotalFeatures: %d\n", TotalFeatures);
  /* for(i = 0; i < count; i++)
     printf("%d  ", flag[i]);
  
     for(i = 0; i < TotalFeatures; i++){
     for(j = 0;j < DIM; j++)
     printf("%f ", features[i]->array[j]);
     printf("\n");
     }
  */
  fclose(feat_file);
  fclose(scp_file);  
  //Uniform initialization all Gaussian mixtures 
  InitializeGMMs(features, DIM, TotalFeatures, numStates);
}

/******************************************************************************
   InitializeGMMs : Uniform Initialization of all GMMs
   inputs : pointer to all feature vectors(VECTORS_OF_F_VECTORS), dimension of feature vectors, 
   totalNumFeatures, numStates
   outputs : Initialized GMMs with means vector and variance vector, posterior probabilities 
   all feature vectors
******************************************************************************/
void InitializeGMMs(VECTOR_OF_F_VECTORS *features, int Dim, int totalNumFeatures, int *numStates){
  VECTOR_OF_F_VECTORS                  *mixtureMeans, *mixtureVars, *allMixtureMeans, *allMixtureVars;
  VECTOR_OF_F_VECTORS                  *featuresForClustering;
  int                                  *numMixEachState;
  float                                *posterior[*numStates];
  int                                  i, j;
  int                                  *mixtureElemCount[*numStates];
  int                                  VQIter = 10, GMMIter = 30;//Number of iterations
  float                                probScaleFactor = 1.0;
  int                                  ditherMean = 1;
  int                                  varianceNormalize = 1;
  int                                  *numElemEachState;
  float                                *Pi;
  //allocate memory  
  featuresForClustering = (VECTOR_OF_F_VECTORS *)calloc(totalNumFeatures, sizeof(VECTOR_OF_F_VECTORS));  
  for(i = 0; i < totalNumFeatures; i++){
    featuresForClustering[i] = (F_VECTOR *)AllocFVector(DIM);
  }
  numElemEachState           = (int *)calloc(*numStates, sizeof(int));
  Pi                         = (float *)calloc(*numStates, sizeof(float));
  allMixtureMeans            = (VECTOR_OF_F_VECTORS *)calloc(INITIAL_NUM_MIX * (*numStates), sizeof(VECTOR_OF_F_VECTORS));
  mixtureMeans               = (VECTOR_OF_F_VECTORS *)calloc(INITIAL_NUM_MIX, sizeof(VECTOR_OF_F_VECTORS));
  allMixtureVars             = (VECTOR_OF_F_VECTORS *)calloc(INITIAL_NUM_MIX * (*numStates), sizeof(VECTOR_OF_F_VECTORS));
  mixtureVars                = (VECTOR_OF_F_VECTORS *)calloc(INITIAL_NUM_MIX, sizeof(VECTOR_OF_F_VECTORS));
  numMixEachState            = (int *)calloc((*numStates), sizeof(int));
  
  for(i = 0; i < INITIAL_NUM_MIX * (*numStates); i++){
    allMixtureMeans[i]               = (F_VECTOR *)AllocFVector(Dim);
    allMixtureVars[i]                = (F_VECTOR *)AllocFVector(Dim);
    allMixtureMeans[i]->numElements  = DIM;
    allMixtureVars[i]->numElements   = DIM;
  }
  for(i = 0; i < INITIAL_NUM_MIX * (*numStates); i++){
    mixtureMeans[i]                  = (F_VECTOR *)AllocFVector(Dim);
    mixtureMeans[i]->numElements     = DIM;
    mixtureVars[i]                   = (F_VECTOR *)AllocFVector(Dim);
    mixtureVars[i]->numElements      = DIM;
  }
  for(i = 0; i < (*numStates); i++){
    mixtureElemCount[i] = (int *)calloc(INITIAL_NUM_MIX, sizeof(int));
  }
  
  for(i = 0; i < (*numStates); i++){
    posterior[i] = (float *)calloc(totalNumFeatures, sizeof(float));
  }
  for(i = 0; i < (*numStates); i++)
    for(j = 0; j < totalNumFeatures; j++)
      posterior[i][j] = 0.0;
  
  //allocate initial num of mixtures to each state
  for(i = 0; i < (*numStates); i++){
    numMixEachState[i] = INITIAL_NUM_MIX;
  }
  for(i = 0; i < totalNumFeatures; i++)
    features[i]->numElements = DIM;
  
  // Build GMM for each state  
  for(i = 0; i < (*numStates); i++){    
    printf("performing uniform Initialization of GMM: %d....\n", i);    
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
    ComputeGMM(featuresForClustering,            numFeatures,              mixtureMeans,             mixtureVars,                       mixtureElemCount[i],                INITIAL_NUM_MIX,          VQIter,           GMMIter,                       probScaleFactor,     ditherMean,                         varianceNormalize,        time(NULL));
    int mixCount = 0, k=0;
    //store current mean and variance into all means and variance
    for(j = 0; j < i; j++)
      mixCount += numMixEachState[j];
    for(j = mixCount; j < mixCount + numMixEachState[i]; j++){
      for(k = 0; k < Dim; k++){
	allMixtureMeans[j]->array[k] = mixtureMeans[j-mixCount]->array[k];
	allMixtureVars[j]->array[k] = mixtureVars[j-mixCount]->array[k];
      }
    }
    printf("Initialization complete\n");
  }
  ComputePosteriorProb(features,    posterior,      allMixtureMeans,       allMixtureVars,     numStates,           numMixEachState,  totalNumFeatures);
  //print all mixtureMeans and mixtureVariances
  /*  for(i = 0; i < (*numStates); i++){
      printf("\nState: %d\n", i);
      int mixCount = 0, k=0;
      for(j = 0; j < i; j++)
      mixCount += numMixEachState[j];
      for(j = mixCount; j < mixCount + numMixEachState[i]; j++){
      printf("\nMean:\n");
      for(k = 0; k < Dim; k++){
      printf("%f  ", allMixtureMeans[j]->array[k]);
      }
      printf("\nVar:\n");
      for(k = 0; k < Dim; k++){
      printf("%f  ", allMixtureVars[j]->array[k]);
      }
      }
      printf("\n");
      }
  */
  /* for(i = 0; i < (*numStates); i++){
     printf("Elements in state : %d\n", i);
     int mixCount = numMixEachState[i];
     for(j = 0; j <  mixCount; j++){
     printf("mixture: %d  contains: %d \n", j, mixtureElemCount[i][j]);
     }
     }
  */
  //  ClusteringAndMerging( features,     allMixtureMeans,       allMixtureVars,    numStates,      numMixEachState,       totalNumFeatures,   posterior,       mixtureElemCount,   numElemEachState,    Pi,   featuresForClustering);
}

/******************************************************************************
   ComputePosteriorProb() : compute posterior probabilities
   features - pointer to posterior matrix , pointer to all feature vectors, allMixtureMeans,
   allMixtureVars,   numstates, numMixEachState, totalNumFeatures, 
   posterior - posterior probability matrix, mixtureElemCount - Element count in each mixture

   outputs : compute posterior probability for each feature vector 
******************************************************************************/
void ComputePosteriorProb( VECTOR_OF_F_VECTORS *features, float **posterior, VECTOR_OF_F_VECTORS *allMixtureMeans, VECTOR_OF_F_VECTORS *allMixtureVars, int *numStates, int *numMixEachState, int totalNumFeatures){
  F_VECTOR *mean              = (F_VECTOR *)AllocFVector(DIM);
  F_VECTOR *var               = (F_VECTOR *)AllocFVector(DIM);
  F_VECTOR *x                 = (F_VECTOR *)AllocFVector(DIM);
  int                         k  = 0, mixCount = 0, i = 0, j = 0;
  //compute posterior probabilities
  // compute numElemEachState -- no of elements in each state
  for(i = 0; i < 1/*totalNumFeatures*/; i++){
    for(k = 0; k < DIM; k++)
      x->array[k] = features[i]->array[k];    
    printf("for feature vector: %d\n", i);
    for(j = 0;j < (*numStates); j++){
      k = 0, mixCount=0;
      for(k = 0; k < j; k++)
	mixCount += numMixEachState[j];
      for(k = mixCount; k < mixCount + numMixEachState[i]; k++){
	int d = 0;
	for(d = 0; d < DIM; d++){
	  mean->array[d] = allMixtureMeans[k]->array[d];
	  mean->numElements = DIM;
	  var->array[d] = allMixtureVars[k]->array[d];
	  var->numElements = DIM;
	  x->numElements = DIM;
	}
	printf("\n\n");
	for(d = 0; d < DIM; d++)
	  printf("%f  ", mean->array[d]);
	printf("\n");
	for(d = 0; d < DIM; d++)
	  printf("%f  ", var->array[d]);
	printf("\n");
	for(d = 0; d < DIM; d++)
	  printf("%f  ", x->array[d]);
	printf("\n\n");
	//printf("num: %d\t", x->numElements);
	volatile float priorProb = 2.0;
	float p = 0, probScaleFactor = 1.0;
	p = ComputeProbability(mean, var,priorProb , x, probScaleFactor);
	printf("prob:  %f\t", p);
	posterior[j][i] += p;
      }
      //printf("%f  ", posterior[j][i]);
    }
    printf("\n");
  }
}

/******************************************************************************
   mergingAndClustering() : Clustering of GMMs and Realignment; Merging of states 
   inputs : 
   features - pointer to all feature vectors, allMixtureMeans, allMixtureVars
   numstates, numMixEachState, totalNumFeatures, 
   posterior - posterior probability matrix, mixtureElemCount - Element count in each mixture

   outputs : Perform viterbi realignment, clustering and merging of states or GMMs
******************************************************************************/
void ClusteringAndMerging(VECTOR_OF_F_VECTORS *features,          VECTOR_OF_F_VECTORS *allMixtureMeans,         VECTOR_OF_F_VECTORS *allMixtureVars,               int *numStates,                  int *numMixEachState,                  int totalNumFeatures,         float **posterior,               int **mixtureElemCount,         int *numElemEachState,                 float *Pi,   VECTOR_OF_F_VECTORS *featuresForClustering){
  // local Variable declaration
  int                           VQIter = 10, GMMIter = 20;
  int                           viterbiIter = 5;
  int                           i = 0, j = 0, k = 0;
  float                         *T1[MAX_NUM_STATES];
  int                           *T2[MAX_NUM_STATES];
  int                           *newStateSeq;
  float                         *BIC;
  for(i = 0; i < (*numStates); i++){
    BIC = (float *)calloc((*numStates), sizeof(float));
  }
  for(i = 0; i < (*numStates); i++){
    T1[i] = (float *)calloc(totalNumFeatures, sizeof(float));
    T2[i] = (int *)calloc(totalNumFeatures, sizeof(int));
  }
  // Merge two GMM States or two mixture models
  //perform the steps again till there are no more states to be merged 
  for(i = 0; i < viterbiIter; i++){
    //perform Viterbi Realignment 
    //newStateSeq     = viterbi(features,     numStates,       allMixtureMeans,      allMixtureVars,      numMixEachState,                totalNumFeatures,     posterior,                numElemEachState,        T1,           T2,                  Pi);
    //perform GMM clustering for all states
    int s = 0;
    for(s = 0; s < (*numStates); s++){      
      // find all elements present in state s using viterbi_realignment
      int count = 0, d = 0; // to count number of features in a state
      for(j = 0; j < totalNumFeatures; j++){
	if(newStateSeq[j] == s){
	  for(d = 0; d < DIM; d++){
	    featuresForClustering[count]->array[d] = features[j]->array[d];
	    featuresForClustering[count]->numElements = DIM;
	  }
	  count++;//increment count
	}
      }
      extern VECTOR_OF_F_VECTORS *mixtureMeans, *mixtureVars;
      extern float probScaleFactor; 
      extern int varianceNormalize, ditherMean;
      int numFeatures = count;
      ComputeGMM(featuresForClustering,            numFeatures,              mixtureMeans,             mixtureVars,                       mixtureElemCount[i],                numMixEachState[s],          VQIter,           GMMIter,                  probScaleFactor,       ditherMean,                         varianceNormalize,        time(NULL));
      int mixCount = 0, k=0;
      //store current mean and variance into all means and variance
      for(j = 0; j < i; j++)
	mixCount += numMixEachState[j];
      for(j = mixCount; j < mixCount + numMixEachState[i]; j++){
	for(k = 0; k < DIM; k++){
	  allMixtureMeans[j]->array[k] = mixtureMeans[j-mixCount]->array[k];
	  allMixtureVars[j]->array[k] = mixtureVars[j-mixCount]->array[k];
	}
      }
    }//GMM Clustering for each state has completed 
    //calculate posterior probabilities
    ComputePosteriorProb(features,    posterior,      allMixtureMeans,       allMixtureVars,     numStates,       numMixEachState,      totalNumFeatures);
  }
  // Calculate BIC Value between Each state   
  CalculateBIC(BIC,    numStates,    posterior,      totalNumFeatures,    newStateSeq,    numMixEachState , numElemEachState);
  // find minimum delta BIC states 
  
  // Merge two states 
  
}

/******************************************************************************
   CalculateBIC() : Calculate BIC value between all the states 
   input:    posterior probability matrix, stateSequence matrix , numMixEachState- to estimate free parameters, numStates, 
             
   outputs : compute BIC for each state
******************************************************************************/
void CalculateBIC(float *BIC, int *numStates,  float **posterior, int totalNumFeatures, int *stateSeq, int *numMixEachState,               int *numElemEachState){

  int                       i = 0, j = 0, k = 0;
  for(i = 0; i < *numStates; i++)
    numElemEachState[i] = 0;
  
  for(i = 0; i < totalNumFeatures; i++){
    int s = stateSeq[i];
    BIC[s] += posterior[s][i];
    numElemEachState[s] += 1;
  }
  for(i = 0; i < *numStates; i++){
    float penalty = 0;
    if(numElemEachState[i] > 0){      
      penalty = numMixEachState[i]  * log(numElemEachState[i]);
    }
    else
      printf("singularity detected: %d state have zero elements\n", i);
    BIC[i] = BIC[i] - penalty;
  }
}
