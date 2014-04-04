#include "../include/viterbi_realign.h"
/******************************************************************************
   viterbi_realign : Perform Viterbi Realignment on given observation vectors
   inputs : Observation Space (pointer to all training feature vectors) same as Sequence of Observations, 
   allMixtureMeans of GMMs - pointer to means of gmms, 
   allMixtureVars of GMMs, (*numStates), numMixEachState, Transition matrix (not using currently), 
   Emission probability matrix Bj(Ot), Initial probabilities (Pi)

   outputs : Most likely hidden state sequence to observation for each input observation vector
******************************************************************************/

int* viterbi( VECTOR_OF_F_VECTORS *features,           int *numStates,                    VECTOR_OF_F_VECTORS *allMixtureMeans, 
	      VECTOR_OF_F_VECTORS *allMixtureVars,     int totalNumFeatures,              float **B,        int *numElemEachState,
	      float *T1[],                             int *T2[],                         float *Pi, int *hiddenStateSeq){
  
  //use numElemEachState for minimum duration constraint
  int                              i = 0, j = 0, k = 0;
  float A                          = (-1) * log(*numStates); // uniform transition probabilities
  int                              max_prob_state;
  //clear any remainings
  for(i = 0; i < *numStates; i++){
    for(j = 0; j < totalNumFeatures; j++){
      T1[i][j] = 0.0; //T1 as delta matrix(in book)
      T2[i][j] = 0;  // T2 as Si matrix
    }
  }
  for(i = 0; i < (*numStates); i++){
    //printf("Pi[%d]: %f\n", i, Pi[i]);
    T1[i][0] = Pi[i] + B[i][0];
    T2[i][0] = 0;
    //printf("T1[%d][0]: %f\n", i, T1[i][0]);
  }
  printf("A: %f \n", A);
  for(i = 1; i < totalNumFeatures; i++){
    //printf("\n%d:  ", i);
    for(j = 0; j < (*numStates); j++){
      float max_prob = -9999999.0;
      int max_prob_state = 0;
      //printf(" %f  ", B[j][i]);
      for(k = 0; k < (*numStates); k++){
	if((T1[k][i-1] + A ) > max_prob){
	  max_prob = T1[k][i-1] + A;
	  max_prob_state = k;
	}
      }
      T1[j][i] = max_prob + B[j][i];
      T2[j][i] = max_prob_state;
      //printf("max_prob_state: %d\n", max_prob_state);
    }
  }
  int *Z = (int *)calloc(totalNumFeatures, sizeof(int));
  // there is an logical algorithmic error in below code 
  float max_prob = -999999;
  int max_state = 0;
  for(i = 0; i < *numStates; i++){
    if( T1[i][totalNumFeatures-1] > max_prob){
      max_prob = T1[i][totalNumFeatures-1];
      max_state = i;
    }
  }
  hiddenStateSeq[totalNumFeatures-1] = max_state;
  int t = 0;
  for(t = totalNumFeatures-2; t >= 0; t--){
    int idx = hiddenStateSeq[t+1];
    hiddenStateSeq[t] = T2[idx][t+1];
  }
  return hiddenStateSeq;
  /*return CheckMinDuration(features , hiddenStateSeq,           numStates,               allMixtureMeans,            allMixtureVars, 
			  totalNumFeatures,        B,                          numElemEachState,
			  T1,                       T2,                      Pi);*/
}


/******************************************************************************
   Checkforminimumduration() : 
   inputs : Observation Space (pointer to all training feature vectors) same as Sequence of Observations, 
   allMixtureMeans of GMMs - pointer to means of gmms, 
   allMixtureVars of GMMs, (*numStates), numMixEachState, Transition matrix (not using currently), 
   Emission probability matrix Bj(Ot), Initial probabilities (Pi)

   outputs : check for minimum duration and drops any state which has less than MIN_DUR frames 
******************************************************************************/
int* CheckMinDuration(    int *hiddenStateSeq,                   int *numStates, 
			  int numFeatures,		      float **B){
  int                 i = 0, j = 0, s = 0;
  float *segProb = (float *)calloc(*numStates, sizeof(float ));
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
