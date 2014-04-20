/*
Written By: Jitendra Prasad Keer
BTech, CSE, IIT Mandi
 */

#include "../../include/preprocessing.h"
#include "../../include/config_multi.h"
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
  int                                  States = INITIAL_STATES;
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
  i = -1; // count file feature vectors
  while(!feof(feat_file) && !ferror(feat_file)){
    i++;
    float number = 0;
    for(j = 0; j < DIM; j++){
      fscanf(feat_file, "%f", &number);
      feat->array[j] = number;
    }
    //check if i is voiced feature vector or not
    if(flag[i] == 0){
      //unvoiced frame
      featureSpace[i][0] = 0;//counter
      featureSpace[i][1] = -1;//counter in voiced features
      featureSpace[i][2] = -1;//speaker number
      continue;
    }
    else{
      //this feature vector is voiced 
      featureSpace[i][0] = 1;
      featureSpace[i][1] = count;
      for(j = 0; j < DIM; j++)
	features[count]->array[j] = feat->array[j];
      features[count]->numElements = DIM;
      count++;
      //      printf("%d  ", i);
    }
  }
  TotalFeatures = count;
  allFeaturesCount = i;
  printf("TotalFeatures: %d  allFeaturesCount: %d\n", TotalFeatures, allFeaturesCount);
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
  // EXTRACT FILE NAME
  char path[512];
  for(i = 0; argv[1][i] != '\0'; i++){
    path[i] = argv[1][i];
  }
  path[i] = '\0';  
  int len = strlen(path);
  printf("path: %s  l: %d\n", path, len);
  int start = 0, end = 0;
  for(i = len-1; i >= 0; i--){
    if(path[i] == '.'){
      end = i;
    }
    if(path[i] == '/'){
      start = i;
      break;
    }
  }
  printf("start: %d  end: %d\n", start, end);
  for(i = start+1; i < end; i++){
    fileName[i-start-1] = path[i];
  }
  fileName[i] = '\0';
  printf("file name: %s \n", fileName);
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
  VECTOR_OF_F_VECTORS                  *allMixtureMeans, *allMixtureVars;

  int                                  *numMixEachState;
  float                                *posterior[*numStates];
  int                                  i, j;
  float                                *mixtureElemCount[*numStates];
  int                                  VQIter = 10, GMMIter = 30;//Number of iterations  
  int                                  *numElemEachState;
  float                                *Pi;
  float                                *mixtureWeight;
  //allocate memory  
  featuresForClustering = (VECTOR_OF_F_VECTORS *)calloc(totalNumFeatures, sizeof(VECTOR_OF_F_VECTORS));
  for(i = 0; i < totalNumFeatures; i++){
    featuresForClustering[i] = (F_VECTOR *)AllocFVector(DIM);
  }
  numElemEachState           = (int *)calloc(*numStates, sizeof(int));
  Pi                         = (float *)calloc(*numStates, sizeof(float));
  allMixtureMeans            = (VECTOR_OF_F_VECTORS *)calloc(INITIAL_NUM_MIX * (*numStates), sizeof(VECTOR_OF_F_VECTORS));
  mixtureMeans               = (VECTOR_OF_F_VECTORS *)calloc(MAX_NUM_MIX, sizeof(VECTOR_OF_F_VECTORS));
  allMixtureVars             = (VECTOR_OF_F_VECTORS *)calloc(INITIAL_NUM_MIX * (*numStates), sizeof(VECTOR_OF_F_VECTORS));
  mixtureVars                = (VECTOR_OF_F_VECTORS *)calloc(MAX_NUM_MIX, sizeof(VECTOR_OF_F_VECTORS));
  numMixEachState            = (int *)calloc((*numStates), sizeof(int));
  mixtureWeight              = (float *)calloc(MAX_NUM_MIX, sizeof(float));
  
  for(i = 0; i < MAX_NUM_MIX * (*numStates); i++){
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
  int total_mix = (*numStates) * INITIAL_NUM_MIX;
  for(i = 0; i < (*numStates); i++){
    mixtureElemCount[i] = (float *)calloc(total_mix, sizeof(float));
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
    ComputeGMM( featuresForClustering,            numFeatures,              mixtureMeans,             mixtureVars,
		mixtureElemCount[i],              INITIAL_NUM_MIX,          VQIter,                   GMMIter,             
		probScaleFactor,                  ditherMean,               varianceNormalize,        time(NULL));
    int mixCount = 0, k=0;
    //store current mean and variance into all means and variance
    for(j = 0; j < i; j++)
      mixCount += numMixEachState[j];
    for(j = mixCount; j < mixCount + numMixEachState[i]; j++){
      for(k = 0; k < Dim; k++){
	allMixtureMeans[j]->array[k] = mixtureMeans[j-mixCount]->array[k];
	allMixtureVars[j]->array[k] = mixtureVars[j-mixCount]->array[k];
      }
      float mix_wt = mixtureElemCount[i][j-mixCount] / numElemEachState[i];
      mixtureWeight[j] = mix_wt;
    }
    printf("Initialization complete\n\n");
  }
  ComputePosteriorProb(features,    posterior,      allMixtureMeans,       allMixtureVars,     numStates,    
		       numMixEachState,  totalNumFeatures,                 mixtureWeight);
  // INitial state probabilities 
  for(i = 0; i < (*numStates); i++){
    Pi[i] = log((float)1.0/(*numStates));
    //printf("Pi[%d]: %f\n", i, Pi[i]);
  }
  PrintAllDetails(numStates, numMixEachState, numElemEachState, 
		  mixtureElemCount, allMixtureMeans, mixtureWeight);
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
    }*/
  ClusteringAndMerging( features,     allMixtureMeans,       allMixtureVars,    numStates,      numMixEachState, 
  			totalNumFeatures,   posterior,       mixtureElemCount,   numElemEachState,    Pi,
			mixtureWeight);
}


/******************************************************************************
   mergingAndClustering() : Clustering of GMMs and Realignment; Merging of states 
   inputs : 
   features - pointer to all feature vectors, allMixtureMeans, allMixtureVars
   numstates, numMixEachState, totalNumFeatures, 
   posterior - posterior probability matrix, mixtureElemCount - Element count in each mixture

   outputs : Perform viterbi realignment, clustering and merging of states or GMMs
******************************************************************************/
void ClusteringAndMerging(VECTOR_OF_F_VECTORS *features,                     VECTOR_OF_F_VECTORS *allMixtureMeans,         
			  VECTOR_OF_F_VECTORS *allMixtureVars,               int *numStates,                  int *numMixEachState,
			  int totalNumFeatures,                              float **posterior,          float **mixtureElemCount,
			  int *numElemEachState,                 float *Pi,  float *mixtureWeight)
{
  // local Variable declaration
  int                           VQIter = 7, GMMIter = 23, mergeIter=0, maxMergeIter = 16;
  int                           viterbiIter = 5;
  int                           i = 0, j = 0, k = 0,s= 0,  flag = 1;
  float                         *T1[MAX_NUM_STATES];
  int                           *T2[MAX_NUM_STATES];
  int                           *newStateSeq;
  float                         *deltaBIC[MAX_NUM_STATES];
  for(i = 0; i < (*numStates); i++){
    deltaBIC[i] = (float *)calloc((*numStates), sizeof(float));
  }
  for(i = 0; i < (*numStates); i++){
    T1[i] = (float *)calloc(totalNumFeatures, sizeof(float));
    T2[i] = (int *)calloc(totalNumFeatures, sizeof(int));
  }
  // Merge two GMM States or two mixture models
  //perform the steps repeatedly till there are no more states to be merged 
  while(flag && mergeIter < maxMergeIter){
    printf("\nClustering and Merging Iteration: %d...........\n", mergeIter);
    printf("starting with, States: %d \n", *numStates);
    mergeIter++;
    for(s = 0; s < *numStates; s++)
      printf("state[%d] : %d \n", s, numElemEachState[s]);
    for(i = 0; i < viterbiIter; i++){
      //perform Viterbi Realignment       
      printf("performing viterbi realignment....\n");
      
      newStateSeq     = viterbi(  features , numStates,       allMixtureMeans,      allMixtureVars,      numMixEachState,
				  totalNumFeatures,     posterior,                numElemEachState,        T1,           T2,
				  Pi, mixtureWeight);
      
      printf("Viterbi realignment has successfully completed...\n");
      //find number of elements in each state
      //FindNumberOfElemInEachState(newStateSeq,   numStates, totalNumFeatures,  numElemEachState, Pi);
      
      //perform GMM clustering for all states
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
	//printf("count: %d  num: %d\n", count, numElemEachState[s]);
	int numMix = numMixEachState[s];
	// Cluster using GMM
	ComputeGMM(featuresForClustering,            numFeatures,              mixtureMeans,             mixtureVars,
		   mixtureElemCount[s],     numMix   ,          VQIter,           GMMIter,               probScaleFactor,
		   ditherMean,         varianceNormalize,        time(NULL));
	
	int mixCount = 0, k=0;
	//store current mean and variance into all means and variance
	for(j = 0; j < s; j++)
	  mixCount += numMixEachState[j];
	for(j = mixCount; j < mixCount + numMixEachState[s]; j++){
	  for(k = 0; k < DIM; k++){
	    allMixtureMeans[j]->array[k] = mixtureMeans[j-mixCount]->array[k];
	    allMixtureVars[j]->array[k] = mixtureVars[j-mixCount]->array[k];
	  }
	  float mix_wt = mixtureElemCount[s][j-mixCount] / numElemEachState[s];
	  mixtureWeight[j] = mix_wt;
	}
      }//GMM Clustering for each state has completed 
      //calculate posterior probabilities      
      PrintAllDetails(numStates, numMixEachState, numElemEachState, 
		      mixtureElemCount, allMixtureMeans, mixtureWeight);
      ComputePosteriorProb(features,    posterior,      allMixtureMeans,       allMixtureVars,     numStates,  
			   numMixEachState,      totalNumFeatures,             mixtureWeight);
    }
    //    PrintAllDetails(numStates, numMixEachState, numElemEachState, 
    //		    mixtureElemCount, allMixtureMeans, mixtureWeight);
    // Calculate BIC Value between Each state   
    CalculateBIC(deltaBIC,    features ,        allMixtureMeans,    allMixtureVars,    numStates,        totalNumFeatures,    
		 newStateSeq,     		 numMixEachState , numElemEachState,   posterior);
    int min_i =0, min_j = 0;
    float min = 99999999;
    // find minimum delta BIC states
    for(j = 0; j < *numStates; j++){
      for(k = j+1; k < *numStates; k++){
	if(deltaBIC[j][k] < min){
	  min = deltaBIC[j][k];
	  min_i = j;
	  min_j = k;
	}
      }
    }
    printf("min_i: %d   min_j: %d   minDBIC: %f\n", min_i, min_j, deltaBIC[min_i][min_j]);
    // if no two states can be merged set the flag to False to terminate process
    if(deltaBIC[min_i][min_j] >= 0){
      printf("We need to stop now ..... Final no of states: %d\n", *numStates);
      flag = 0;
    }
    // Merge two states 
    if(flag)
      MergeTwoStates( min_i,   min_j,   allMixtureMeans,    allMixtureVars,      numStates,     numMixEachState,   
		      numElemEachState,                     Pi,                  totalNumFeatures,
		      features,                  newStateSeq,                    mixtureWeight);
    /*printf("Means after Merging of two states.................\n");
      PrintAllDetails(numStates, numMixEachState, numElemEachState, 
      mixtureElemCount, allMixtureMeans);*/
    for(i = 0; i < *numStates; i++)
      printf("mix[%d]: %d\n", i, numMixEachState[i]);
  }
  //WRITE PLOT_File to plot distribution of each data point
  writePlotFile(posterior, totalNumFeatures, numStates);
  //WRITE RTTM FILE
  writeRTTMFile(newStateSeq, numStates, totalNumFeatures, numElemEachState);
}

/******************************************************************************/
void PrintAllDetails(int *numStates, int *numMixEachState, int *numElemEachState, 
		     float **mixtureElemCount, VECTOR_OF_F_VECTORS *allMixtureMeans, float *mixtureWeight){
  int                i = 0, j = 0, s = 0, d =0;
  for(s = 0; s < *numStates; s++){
    printf("Elems in state: %d       %d\n", s, numElemEachState[s]);
    printf("no of mixtures in state: %d     %d \n", s, numMixEachState[s]);
    for(i = 0; i < numMixEachState[s]; i++)
      printf("\tmixture %d : %f\n", i, mixtureElemCount[s][i]);
    //print means of all mixtures
    int mixCount = 0;
    for(j = 0; j < s; j++)
      mixCount += numMixEachState[j];
    for(j = mixCount; j < mixCount + numMixEachState[s]; j++){
      //for(d = 0; d < DIM; d++)
      //printf("%f ", allMixtureMeans[j]->array[d]);
      printf("mix:%d  %f\n", j, mixtureWeight[j]);
      printf("\n");
    }
  }
}
/******************************************************************************/
void writePlotFile(float **posterior, int totalNumFeatures, int *numStates){
  int                    i = 0, j = 0, s= 0, d = 0;
  int                    *featureState;
  FILE                   *plotFile;
  featureState = (int *)calloc(totalNumFeatures, sizeof(int *));
  for(i = 0; i < totalNumFeatures; i++){
    float max = -999999;
    int max_s = 0;
    for(s = 0; s < *numStates; s++){
      if(posterior[s][i] > max){
	max = posterior[s][i];
	max_s = s;
      }
    }
    featureState[i] = max_s;
  }
  plotFile = fopen("plot_data.txt", "w");
  if(!plotFile){
    fprintf(stderr, "unable to open plot file\n");
  }
  for(i = 0; i < totalNumFeatures; i++){
    fprintf(plotFile, "%d    %d\n", i, featureState[i]);
  }
  fclose(plotFile);
}

/******************************************************************************
   writeRTTMFile() : write rttm file in format given by NIST
   input:    state Sequence, no of states, total no of features, no of elements in each state
             
   outputs : rtt file
******************************************************************************/
void writeRTTMFile(int *stateSeq, int *numStates, int totalNumFeatures, int *numElemEachState){
  FILE                            *rttm;
  int                             i = 0, j = 0, k = 0, s = 0;
  //open file to write
  char file[100];
  strcpy(file, fileName);
  rttm = fopen(strcat(file, ".rttm"), "w");
  if(!rttm){
    fprintf(stderr, "Unable to open rttm file for writing\n");
    exit(0);
  }
  for(s = 0; s < *numStates; s++){
    char c = s + '0';
    char spkr[4] = {'_', c , '_', '\0'};
    char spkr_name[128];
    strcpy(spkr_name, fileName);
    strcat(spkr_name, spkr);
    printf("spkr_name: %s\n", spkr_name);
    fprintf(rttm, "SPKR-INFO %s 1 <NA> <NA> <NA> unknown %s <NA>\n", fileName, spkr_name);
    int start = 0, dur = 0;
    for(k = 0; k < allFeaturesCount; k++){
      int index = 0;
      index = featureSpace[k][1];
      //if frame is voiced
      if(featureSpace[k][0] && index >= 0){
	if(stateSeq[index] != s && dur > 0){
	  //write into the file
	  fprintf(rttm, "SPEAKER %s 1 %f %f <NA> <NA> %s <NA>\n", fileName, (float)start * 0.010, (float) dur * 0.010, spkr_name);
	  dur = 0;
	}
	if(stateSeq[index] == s && dur == 0){
	  start = k;
	}
	if(stateSeq[index] == s){
	  dur++;
	}
      }
    }
  }
  //done

}
/******************************************************************************
   MergeTwoStates() : Merge two states 
   input:    min_i_state , min_j_state , allMeans  , all Vars , number of states, no of mix in each state
             
   outputs : Two states merged and reduction in number of states by one
******************************************************************************/
void MergeTwoStates( int s_i,            int s_j,              VECTOR_OF_F_VECTORS *allMixtureMeans,     
		     VECTOR_OF_F_VECTORS *allMixtureVars,      int *numStates,                        int *numMixEachState,  
		     int *numElemEachState,                    float *Pi,                             int totalNumFeatures,
		     VECTOR_OF_F_VECTORS *features,            int *stateSeq,                         float *mixtureWeight)
{
  printf("merging state %d and %d........\n", s_i, s_j);
  int                                  i = 0, j = 0, d = 0, s = 0, min_s = 0, max_s=0;
  VECTOR_OF_F_VECTORS                  *tempAllMixtureMeans, *tempAllMixtureVars;
  tempAllMixtureMeans   = (VECTOR_OF_F_VECTORS *) calloc(MAX_NUM_MIX * (*numStates) , sizeof(VECTOR_OF_F_VECTORS));
  tempAllMixtureVars    = (VECTOR_OF_F_VECTORS *) calloc(MAX_NUM_MIX * (*numStates) , sizeof(VECTOR_OF_F_VECTORS));
  for(i = 0; i < MAX_NUM_MIX * (*numStates); i++){
    tempAllMixtureMeans[i] = (F_VECTOR *) AllocFVector(DIM);
    tempAllMixtureVars[i]  = (F_VECTOR *) AllocFVector(DIM);
  }
  //find smaller state
  if(s_i < s_j){
    min_s = s_i;
    max_s = s_j;
  }
  else{
    min_s = s_j;
    max_s = s_i;
  }
  // COPY ORIGINAL MEANS AND VARS INTO TEMPS
  int totalMix = 0;
  for(s = 0; s < *numStates; s++)
    totalMix += numMixEachState[s];
  for(i = 0; i < totalMix; i++){
    for(d = 0; d < DIM; d++){
      tempAllMixtureMeans[i]->array[d] = allMixtureMeans[i]->array[d];
      tempAllMixtureMeans[i]->numElements = DIM;
      tempAllMixtureVars[i]->array[d]  = allMixtureVars[i]->array[d];
      tempAllMixtureMeans[i]->numElements = DIM;
    }
  }
  //copy all means and vars till min_s
  /*for(s = 0; s < min_s; s++){
    for(d = 0; d < DIM; d++){
    allMixtureMeans[i]->array[d] = tempAllMixtureMeans[i]->array[d];
    allMixtureMeans[i]->numElements = DIM;
    allMixtureVars[i]->array[d]  = tempAllMixtureVars[i]->array[d];
    allMixtureMeans[i]->numElements = DIM;
    }
    }*/
  // copy max_s mixture means and vars also into min_s state  
  int mixCount = 0;
  for(s = 0; s <= min_s; s++)
    mixCount += numMixEachState[s];
  int go_to_max = 0;
  for(s = 0; s < max_s; s++)
    go_to_max += numMixEachState[s];
  
  for(i = mixCount, j = go_to_max; i < mixCount + numMixEachState[max_s], 
	j < go_to_max + numMixEachState[max_s]; i++, j++){
    for(d = 0; d < DIM; d++){
      allMixtureMeans[i]->array[d] = tempAllMixtureMeans[j]->array[d];
      allMixtureMeans[i]->numElements = DIM;
      allMixtureVars[i]->array[d]  = tempAllMixtureVars[j]->array[d];
      allMixtureMeans[i]->numElements = DIM;
    }
  }  
  numElemEachState[min_s] += numElemEachState[max_s];
  // decrease number of states by one 
  *numStates = *numStates - 1;
  for(s = min_s+1; s < *numStates - 1; s++){
    if(s < max_s){
      mixCount = 0;
      for(i = 0; i < s; i++)
	mixCount += numMixEachState[i];
      go_to_max = mixCount;
      mixCount += numMixEachState[max_s];      
      
      for(i = mixCount, j = go_to_max; i < mixCount + numMixEachState[s], 
	    j < go_to_max + numMixEachState[s]; i++, j++){
	for(d = 0; d < DIM; d++){
	  allMixtureMeans[i]->array[d] = tempAllMixtureMeans[j]->array[d];
	  allMixtureMeans[i]->numElements = DIM;
	  allMixtureVars[i]->array[d]  = tempAllMixtureVars[j]->array[d];
	  allMixtureMeans[i]->numElements = DIM;
	}
      }
    }
    else if(s >= max_s){
      mixCount = 0;
      for(i = 0; i < s; i++)
	mixCount += numMixEachState[i];
      for(i = mixCount; i < mixCount + numMixEachState[s+1]; i++){
	for(d = 0; d < DIM; d++){
	  allMixtureMeans[i]->array[d] = tempAllMixtureMeans[i]->array[d];
	  allMixtureMeans[i]->numElements = DIM;
	  allMixtureVars[i]->array[d]  = tempAllMixtureVars[i]->array[d];
	  allMixtureMeans[i]->numElements = DIM;
	}
      }
      numMixEachState[s] = numMixEachState[s+1];
      numElemEachState[s] = numElemEachState[s+1];
    }
  }
  numMixEachState[min_s] += numMixEachState[max_s];
  for(s = 0; s < *numStates; s++)
    Pi[s] = (float)log((float) numElemEachState[s]/totalNumFeatures);
  //MODEL MIN_S STATE AGAIN USING FEATURES OF BOTH STATES 
  int count = 0, k = 0;
  for(k = 0; k < totalNumFeatures; k++){
    if( stateSeq[k] == i || stateSeq[k] == j ){
      for( d = 0; d < DIM; d++){
	featuresForClustering[count]->array[d] = features[k]->array[d];
	featuresForClustering[count]->numElements = DIM;
      }
      count++;
    }
  }
  //compute GMM for merged state min_s
  extern VECTOR_OF_F_VECTORS *mixtureMeans, *mixtureVars;
  extern float probScaleFactor; 
  extern int varianceNormalize, ditherMean;
  int numFeatures = count;
  //printf("count: %d  num: %d\n", count, numElemEachState[s]);
  int numMix = numMixEachState[min_s];
  float *mixtureElemCount = (float *)calloc(MAX_NUM_MIX, sizeof(float *));
  // Cluster using GMM
  ComputeGMM(featuresForClustering,            numFeatures,              mixtureMeans,             mixtureVars,
	     mixtureElemCount,     numMix   ,          10,           20,               probScaleFactor,
	     ditherMean,         varianceNormalize,        time(NULL));	
  //store current mean and variance into all means and variance
  for(j = 0; j < min_s; j++)
    mixCount += numMixEachState[j];
  for(j = mixCount; j < mixCount + numMixEachState[min_s]; j++){
    for(k = 0; k < DIM; k++){
      allMixtureMeans[j]->array[k] = mixtureMeans[j-mixCount]->array[k];
      allMixtureVars[j]->array[k] = mixtureVars[j-mixCount]->array[k];
    }
    float mix_wt = mixtureElemCount[j-mixCount] / numElemEachState[min_s];
    mixtureWeight[j] = mix_wt;
  }
  // FREE TEMP VECTORS
  for(i = 0; i < MAX_NUM_MIX * (*numStates); i++){
    free(tempAllMixtureMeans[i]);
    free(tempAllMixtureVars[i]);
  }
  free(tempAllMixtureMeans);
  free(tempAllMixtureVars);
}

/******************************************************************************
   Findnumberofelementsineachstate:Find Number of elements in each state
   input:    newStateSequence , numStates , numElemEachState , Pi 
             
   outputs : calculate number of elements in each state and Pi( Initital state probabilities)
******************************************************************************/
void FindNumberOfElemInEachState(int *newStateSeq, int *numStates, int totalNumFeatures,  int *numElemEachState, float *Pi){
  int                      i = 0, j = 0, k = 0;
  
  printf("\nCalculating number of elements in each state.....\n");
  for(i = 0; i < *numStates; i++)
    numElemEachState[i] = 0;
  for(i = 0; i < totalNumFeatures; i++){
    int s = newStateSeq[i];    
    numElemEachState[s] += 1;
  }
  for(i = 0; i < *numStates; i++){
    printf("Elements in state : %d    %d \n", i, numElemEachState[i]);
  }
  // modify Pi (initial state probabilities)
  for(i = 0; i < *numStates; i++){
    Pi[i] = numElemEachState[i]/totalNumFeatures;
    if(Pi[i] != 0){
      Pi[i] = log(Pi[i]);
    }
  }
}

/******************************************************************************
   ComputePosteriorProb() : compute posterior probabilities
   features - pointer to posterior matrix , pointer to all feature vectors, allMixtureMeans,
   allMixtureVars,   numstates, numMixEachState, totalNumFeatures, 
   posterior - posterior probability matrix, mixtureElemCount - Element count in each mixture

   outputs : compute posterior probability for each feature vector 
******************************************************************************/
void ComputePosteriorProb( VECTOR_OF_F_VECTORS *features, float **posterior, VECTOR_OF_F_VECTORS *allMixtureMeans, 
			   VECTOR_OF_F_VECTORS *allMixtureVars, int *numStates, int *numMixEachState, int totalNumFeatures, 
			   float               *mixtureWeight){
  mean              = (F_VECTOR *)AllocFVector(DIM);
  var               = (F_VECTOR *)AllocFVector(DIM);
  x                 = (F_VECTOR *)AllocFVector(DIM);
  int               k  = 0, mixCount = 0, i = 0, j = 0;
  //compute posterior probabilities
  // compute numElemEachState -- no of elements in each state
  printf("calculating posterior probabilities for all feature vectors.......\n");
  for(i = 0; i < totalNumFeatures; i++){
    for(j = 0; j < *numStates; j++){
      posterior[j][i] = 0.0;
    }
  }
  //for(i = 0; i < INITIAL_NUM_MIX * (*numStates); i++)
  //printf("mix %d : %f\n", i, mixtureWeight[i]);
  for(i = 0; i < totalNumFeatures; i++){
    for(k = 0; k < DIM; k++)
      x->array[k] = features[i]->array[k];
    //printf("\nfor feature vector: %d\n", i);
    for(j = 0;j < (*numStates); j++){//for each state
      k = 0, mixCount=0;
      for(k = 0; k < j; k++)
	mixCount += numMixEachState[k];
      for(k = mixCount; k < mixCount + numMixEachState[j]; k++){
	int d = 0;
	for(d = 0; d < DIM; d++){
	  mean->array[d] = allMixtureMeans[k]->array[d];
	  mean->numElements = DIM;
	  var->array[d] = allMixtureVars[k]->array[d];
	  var->numElements = DIM;
	  x->numElements = DIM;
	}
	volatile float priorProb = 1.0;
	float p = 0, probScaleFactor = 1.0;
	p = ComputeProbability(mean, var,priorProb , x, probScaleFactor);
	float wt = mixtureWeight[k];
	//printf("prob:  %f\t", p);
	posterior[j][i] = posterior[j][i] + wt * p;
      }
      //printf("%f  ", posterior[j][i]);
    }
    //printf("\n");
  }
  printf("finished...calculation of posterior probabilities.\n");
}

/******************************************************************************
   BIC_Modified() : Calculate BIC value between all the states 
   input:    posterior probability matrix, stateSequence matrix , numMixEachState- to estimate free parameters, numStates, 
             
   outputs : compute BIC for each state
******************************************************************************/
void CalculateBIC( float **deltaBIC,   VECTOR_OF_F_VECTORS *features,    VECTOR_OF_F_VECTORS *allMixtureMeans, 
		   VECTOR_OF_F_VECTORS *allMixtureVars,                  int *numStates,   int totalNumFeatures, int *stateSeq,  
		   int *numMixEachState,               int *numElemEachState, float **posterior)
{
  int                               s = 0, i = 0, j = 0, k = 0, count = 0, d= 0;
  float                             L0 = 0, L1 = 0;
  //extern VECTOR_OF_F_VECTORS        *mixtureMeans, *mixtureVars, *featuresForClustering;
  
  for(i = 0; i < *numStates; i++){
    for(j = i+1; j < *numStates; j++){
      //calculate Delta BIC between state i and state j       
      L0 = 0, L1 = 0;
      //calculate L1 logSum of probs from original models 
      count = 0; int checker = 0;
      for(k = 0; k < totalNumFeatures; k++){
	if( stateSeq[k] == i || stateSeq[k] == j ){
	  for( d = 0; d < DIM; d++){
	    featuresForClustering[count]->array[d] = features[k]->array[d];
	    featuresForClustering[count]->numElements = DIM;
	  }
	  count++;
	}
	if( stateSeq[k] == i ){
	  L1 += posterior[i][k];
	  if(checker < 100){
	    printf("p: %f\t", posterior[i][k]);
	    checker++;
	  }
	}
	else if( stateSeq[k] == j){
	  L1 += posterior[j][k];
	  if(checker < 100){
	    printf("p: %f\t", posterior[j][k]);
	    checker++;
	  }
	}
      }
      checker = 0;
      //printf("sum of log likelihood from both models L1: %f\n", L1);
      //Calculate L0
      // Model elements of both these models using single gaussian        	  
      extern float probScaleFactor; 
      extern int varianceNormalize, ditherMean;
      int numFeatures = count;
      int numMixtures = numMixEachState[i] + numMixEachState[j];
      float *mixtureElemCount = (float *) AllocFloatArray(mixtureElemCount, 
							  numMixtures);
      int VQIter = 10, GMMIter = 20;
      ComputeGMM(featuresForClustering,            numFeatures,              mixtureMeans,             mixtureVars,
		 mixtureElemCount,                numMixtures,          VQIter,           GMMIter,                 
		 probScaleFactor,       ditherMean,                         varianceNormalize,        time(NULL));
      
      // now calculate likelihood of each feature from newly obtained model
      extern F_VECTOR *mean, *var, *x;
      int mix = 0;
      for(k = 0; k < numFeatures; k++){
	float P = 0;
	for(mix = 0; mix < numMixtures; mix++){
	  for(d = 0; d < DIM; d++){
	    x->array[d] = featuresForClustering[k]->array[d];
	    x->numElements = DIM;
	    mean->array[d] = mixtureMeans[mix]->array[d];
	    mean->numElements = DIM;
	    var->array[d] = mixtureVars[mix]->array[d];
	    var->numElements = DIM;
	  }
	  float p = 0;
	  extern float probScaleFactor;
	  float priorProb = 1.0;
	  p = ComputeProbability(mean, var, priorProb, x, probScaleFactor);
	  float wt = mixtureElemCount[mix] / numFeatures;
	  P = P + wt * p;
	}
	if(checker < 100){
	  printf("A: %f\t", P);
	  checker++;
	}
	L0 += P;
      }
      //printf("total log sum probability from new Model is L0: %f\n", L0);
      // Now we have both L0 and L1
      int n1 = numMixEachState[i];
      int n2 = numMixEachState[j];
      int deltaK = 1;      
      float Penalty = 0.5 * LAMBDA * deltaK * log(numFeatures);
      deltaBIC[i][j] = L1 - L0 - Penalty;
      printf("i:%d j:%d  deltaBIC: %f L1: %f  L0: %f\n", i, j , deltaBIC[i][j], L1, L0);
    }
  }
}


/*printf("\n\n");
  for(d = 0; d < DIM; d++)
  printf("%f  ", mean->array[d]);
  printf("\n");
  for(d = 0; d < DIM; d++)
  printf("%f  ", var->array[d]);
  printf("\n");
  for(d = 0; d < DIM; d++)
  printf("%f  ", x->array[d]);
  printf("\n\n");*/
//printf("num: %d\t", x->numElements);

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
