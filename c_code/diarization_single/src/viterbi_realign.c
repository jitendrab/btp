#include "../include/viterbi_realign.h"
/******************************************************************************
   viterbi_realign : Perform Viterbi Realignment on given observation vectors
   inputs : Observation Space (pointer to all training feature vectors) same as Sequence of Observations, 
   allMixtureMeans of GMMs - pointer to means of gmms, 
   allMixtureVars of GMMs, (*numStates), numMixEachState, Transition matrix (not using currently), 
   Emission probability matrix Bj(Ot), Initial probabilities (Pi)

   outputs : Most likely hidden state sequence to observation for each input observation vector
******************************************************************************/

int* viterbi( VECTOR_OF_F_VECTORS *features, hmm* mdHMM,   int totalNumFeatures,              float **posterior,     int *numElemEachState,
	      float *delta[],                             int *psi[],                         float *Pi,         int *path){
  
  float **B = (float **)calloc(totalNumFeatures, sizeof(float *));
  int                              numStates = *(mdHMM->numStates);
  int                              i = 0, j = 0, k = 0;
  float A                          = (-1) * log(numStates); // uniform transition probabilities
  int                              max_prob_state;  
  
  for(i = 0; i < totalNumFeatures; i++)
    B[i] = (float *)calloc(numStates * MIN_DUR, sizeof(float ));
  
  //clear any remainings
  for(i = 0; i < numStates; i++){
    for(j = 0; j < totalNumFeatures; j++){
      delta[i][j] = 0.0; //T1 as delta matrix(in book)
      psi[i][j] = 0;  // T2 as Si matrix
    }
  }
  //calculate new B, copy probabilities 
  for(i = 0; i < totalNumFeatures; i++){
    for(j = 0; j < numStates; j++){
      for(k = 0; k < MIN_DUR; k++)
	B[i][j*MIN_DUR + k] = posterior[j][i];
    }
  }
  
  int scale = 0;
  // we are working on log scale
  
  
}


/******************************************************************************
   Checkforminimumduration() : 
   inputs : Observation Space (pointer to all training feature vectors) same as Sequence of Observations, 
   allMixtureMeans of GMMs - pointer to means of gmms, 
   allMixtureVars of GMMs, (*numStates), numMixEachState, Transition matrix (not using currently), 
   Emission probability matrix Bj(Ot), Initial probabilities (Pi)

   outputs : check for minimum duration and drops any state which has less than MIN_DUR frames 
******************************************************************************/

/*int* CheckMinDuration(    int *hiddenStateSeq,                   int *numStates, 
			  int numFeatures,		      float **B){
  int                 i = 0, j = 0, s = 0;
  float *segProb = (float *)calloc(*numStates, sizeof(float *));
  //divide all voiced feature space into MIN_DUR segments
  printf("applying min-duration constraint....\n");
  for(j = 0; j < numFeatures/MIN_DUR; j++){
    for(s = 0; s < *numStates; s++)
      segProb[s] = 0.0;
    
    for(s = 0; s < *numStates; s++){
      for(i = j*MIN_DUR; i < numFeatures, i < (j+1) * MIN_DUR; i++){
	// this i will be a single block or segment
	segProb[s] += B[s][i];
      }
    }
    //find max_state 
    double max = -9999999;
    int max_s = 0;
    for(s = 0; s < *numStates; s++){
      if(segProb[s] > max){
	max = segProb[s];
	max_s = s;
      }
    }
    //max_s is the state where this segment should belong
    for(i = j*MIN_DUR; i < numFeatures, i < (j+1) * MIN_DUR; i++){
      // this i will be a single block or segment
      hiddenStateSeq[i] = max_s;
    }
  }
  free(segProb);
  return hiddenStateSeq;
}
*/

int* CheckMinDuration(    int *hiddenStateSeq,                   int *numStates, 
			  int numFeatures,		      float **B){
  int                 i = 0, j = 0, s = 0;
  int *counter = (int *)calloc(*numStates, sizeof(int ));
  
  //divide all voiced feature space into MIN_DUR segments
  printf("applying min-duration constraint....\n");
  for(j = 0; j < (numFeatures/MIN_DUR) +1; j++){
    for(s = 0; s < *numStates; s++)
      counter[s] = 0;
    
    for(s = 0; s < *numStates; s++){
      for(i = j*MIN_DUR; i < numFeatures, i < (j+1) * MIN_DUR; i++){
	// this i will be a single block or segment
	int idx = hiddenStateSeq[i];
	counter[idx]++;
      }
    }
    //assign this cluster to the maximum 
    int max = -99999, max_s = 0;
    for(s = 0; s < *numStates; s++){
      if(counter[s] > max){
	max = counter[s];
	max_s = s;
      }
    }
    for(i = j*MIN_DUR; i < numFeatures, i < (j+1) * MIN_DUR; i++){
      // this i will be a single block or segment
      hiddenStateSeq[i] = max_s;
    }
  }
  return hiddenStateSeq;
}
