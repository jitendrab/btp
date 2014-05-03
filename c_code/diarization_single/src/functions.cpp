/*
Written By: Jitendra Prasad Keer
BTech, CSE, IIT Mandi
 */

/*-------------------------------------------------------------------------
 *  AllocFloatArray -- Allocates an array of floats
 *    Args:	Array, size of array
 *    Returns:	allocated array
 *    Bugs:	
 * -------------------------------------------------------------------------*/

float * AllocFloatArray(float *array, int npts)
{
  array = (float *) calloc (npts, sizeof(float));
  if (array == NULL) {
    printf("unable to allocate Float array \n");
    exit(-1);
  }
  return(array);
}	/*  End of AllocFloatArray */	

/*-------------------------------------------------------------------------
 *  AllocIntArray -- Allocates an array of Ints
 *    Args:	Array, size of array
 *    Returns:	allocated array
 *    Bugs:	
 * -------------------------------------------------------------------------*/

int *AllocIntArray(int *array, int npts)
{
  array = (int *) calloc(npts,sizeof(int));
  if (array == NULL) {
    printf("unable to allocate Int array \n");
    exit(-1);
  }
  return(array);
}	/*  End of AllocIntArray  */


/*--------------------------------------------------------------------------
    Determine if most elements in an F_VECTOR are zero 
    ---------------------------------------------------------------------------*/
int ZeroFVector (F_VECTOR *fvect) {
  int i = 0;
  int flag = 0;
  int cnt = 0;
  int numElem = ceilf((float)fvect->numElements/3.0);
  while ((i < fvect->numElements) && (cnt < numElem)) {
    flag = (fvect->array[i] == 0);
    i++;
    if (flag)
      cnt++;
  }
  flag = cnt >= numElem;
  return(flag);
}

/******************************************************************************
   AllocFVector : Allocate space for a variable of type F_VECTOR with an 
   array size of numPts.
   inputs : npts
   outputs : an F_VECTOR of size npts
******************************************************************************/

F_VECTOR *AllocFVector(int npts) {
  F_VECTOR *fVect;
  fVect = (F_VECTOR *) malloc(1*sizeof(F_VECTOR));
  if (fVect == NULL) {
    printf("unable to allocate FVector \n");
    exit(-1);
  }
  fVect->numElements = npts;
  fVect->array = (float *) calloc(npts, sizeof(float));
  if (fVect->array == NULL) {
    printf("unable to allocate FVector array \n");
    exit(-1);
  } else
    return(fVect);
}

/******************************************************************************
   CalculateBIC() : Calculate BIC value between all the states 
   input:    posterior probability matrix, stateSequence matrix , numMixEachState- to estimate free parameters, numStates, 
             
   outputs : compute BIC for each state
******************************************************************************/
void CalculateBIC(float **deltaBIC,   VECTOR_OF_F_VECTORS *features,    int *numStates,  float **posterior, int totalNumFeatures, int *stateSeq, int *numMixEachState,               int *numElemEachState){
  int                               s = 0, i = 0, j = 0, k = 0, count = 0, d= 0;
  float                             L0 = 0, L1 = 0;
  //extern VECTOR_OF_F_VECTORS        *mixtureMeans, *mixtureVars, *featuresForClustering;
  
  for(i = 0; i < *numStates; i++){
    for(j = i+1; j < *numStates; j++){     
      //calculate Delta BIC between state i and state j       
      L0 = 0, L1 = 0;
      //calculate L1 logSum of probs from original models 
      count = 0;
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
	  //printf("p: %f\t", posterior[i][k]);
	}
	else if( stateSeq[k] == j){
	   L1 += posterior[j][k];
	  //printf("p: %f\t", posterior[j][k]);
	}
      }
      //printf("sum of log likelihood from both models L1: %f\n", L1);
      //Calculate L0
      // Model elements of both these models using single gaussian        	  
      extern float probScaleFactor; 
      extern int varianceNormalize, ditherMean;
      int numFeatures = count;
      int numMixtures = 1;
      float *mixtureElemCount = (float *) AllocFloatArray(mixtureElemCount, 
							  numMixtures);
      int VQIter = 10, GMMIter = 20;
      ComputeGMM(featuresForClustering,            numFeatures,              mixtureMeans,             mixtureVars,                       mixtureElemCount,                numMixtures,          VQIter,           GMMIter,                  probScaleFactor,       ditherMean,                         varianceNormalize,        time(NULL));
      
      // now calculate likelihood of each feature from newly obtained model
      extern F_VECTOR *mean, *var, *x;
      for(d = 0; d < DIM; d++){
	mean->array[d] = mixtureMeans[0]->array[d];
	mean->numElements = DIM;
	var->array[d] = mixtureVars[0]->array[d];
	var->numElements = DIM;
	x->numElements = DIM;
      }
      for(k = 0; k < numFeatures; k++){
	for(d = 0; d < DIM; d++)
	  x->array[d] = featuresForClustering[k]->array[d];    
	volatile float priorProb = 1.0;
	float p = 0; extern float probScaleFactor;
	p = ComputeProbability(mean, var,priorProb , x, probScaleFactor);
	//printf("p: %f\n", p);
	L0 += p;
      } 
      //printf("total log sum probability from new Model is L0: %f\n", L0);
      // Now we have both L0 and L1
      int n1 = numMixEachState[i];
      int n2 = numMixEachState[j];
      int deltaK = 2 * (n1 + n2 ) - 2;      
      float Penalty = 0.5 * LAMBDA * deltaK * log(numFeatures);
      deltaBIC[i][j] = L0 - L1 - Penalty;
      //printf("i: %d   j: %d   deltaBIC: %f\n\n\n", i, j , deltaBIC[i][j]);
    }
  }
}

/******************************************************************************
   MergeTwoStates() : Merge two states 
   input:    min_i_state , min_j_state , allMeans  , all Vars , number of states, no of mix in each state
             
   outputs : Two states merged and reduction in number of states by one
******************************************************************************/
void MergeTwoStates( int s_i,  int s_j  , VECTOR_OF_F_VECTORS *allMixtureMeans,     
		     VECTOR_OF_F_VECTORS *allMixtureVars,      int *numStates,    int *numMixEachState,  
		     int *numElemEachState,                    float *Pi,          int totalNumFeatures){
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
  for(i = mixCount; i < mixCount + numMixEachState[max_s]; i++){
    for(d = 0; d < DIM; d++){
      allMixtureMeans[i]->array[d] = tempAllMixtureMeans[i]->array[d];
      allMixtureMeans[i]->numElements = DIM;
      allMixtureVars[i]->array[d]  = tempAllMixtureVars[i]->array[d];
      allMixtureMeans[i]->numElements = DIM;
    }
  }
  numMixEachState[min_s] += numMixEachState[max_s];
  numElemEachState[min_s] += numElemEachState[max_s];
  for(s = min_s+1; s < *numStates - 1; s++){    
    if(s < max_s){
      mixCount = 0;
      for(i = 0; i < s; i++)
	mixCount += numMixEachState[i];   
      for(i = mixCount; i < mixCount + numMixEachState[s]; i++){
	for(d = 0; d < DIM; d++){
	  allMixtureMeans[i]->array[d] = tempAllMixtureMeans[i]->array[d];
	  allMixtureMeans[i]->numElements = DIM;
	  allMixtureVars[i]->array[d]  = tempAllMixtureVars[i]->array[d];
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
  // decrease number of states by one 
  *numStates = *numStates - 1;
  for(s = 0; s < *numStates; s++)
    Pi[s] = (float)log((float) numElemEachState[s]/totalNumFeatures);
  // FREE TEMP VECTORS
  for(i = 0; i < MAX_NUM_MIX * (*numStates); i++){
    free(tempAllMixtureMeans[i]);
    free(tempAllMixtureVars[i]);
  }
  free(tempAllMixtureMeans);
  free(tempAllMixtureVars);
}
