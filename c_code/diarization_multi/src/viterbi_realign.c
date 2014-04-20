/*
Written By: Jitendra Prasad Keer
BTech, CSE, IIT Mandi
 */

#include "../../include/viterbi_realign.h"
/******************************************************************************
   viterbi_realign : Perform Viterbi Realignment on given observation vectors
   inputs : Observation Space (pointer to all training feature vectors) same as Sequence of Observations, 
   allMixtureMeans of GMMs - pointer to means of gmms, 
   allMixtureVars of GMMs, (*numStates), numMixEachState, Transition matrix (not using currently), 
   Emission probability matrix Bj(Ot), Initial probabilities (Pi)

   outputs : Most likely hidden state sequence to observation for each input observation vector
******************************************************************************/

int* viterbi( VECTOR_OF_F_VECTORS *features, int *numStates, VECTOR_OF_F_VECTORS *allMixtureMeans, 
	      VECTOR_OF_F_VECTORS *allMixtureVars, int *numMixEachState, int totalNumFeatures, float **B, 
	      int *numElemEachState, float *T1[], int *T2[], float *Pi, float *mixtureWeight){
  
  //use numElemEachState for minimum duration constraint
  int                              i = 0, j = 0, k = 0;
  float A                          = (-1) * log(*numStates); // uniform transition probabilities
  int                              *hiddenStateSeq;
  int                              max_prob_state;
  hiddenStateSeq = (int *)calloc(totalNumFeatures, sizeof(int));
  //clear any remainings
  for(i = 0; i < *numStates; i++){
    for(j = 0; j < totalNumFeatures; j++){
      T1[i][j] = 0.0;
      T2[i][j] = 0;
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
	if(( T1[k][i-1] + A + B[j][i]) > max_prob){
	  max_prob = T1[k][i-1] + A + B[j][i];
	  max_prob_state = k;
	}
      }
      T1[j][i] = max_prob;
      T2[j][i] = max_prob_state;
      //printf("max_prob_state: %d\n", max_prob_state);
    }
  }
  int *Z = (int *)calloc(totalNumFeatures, sizeof(int));
  
  for(i = 0; i < totalNumFeatures; i++){
    float max_prob = -9999999;
    max_prob_state = 0;
    for(k = 0; k < (*numStates); k++){
      if( T1[k][i] > max_prob){
	max_prob = T1[k][i];
	max_prob_state = k;
      }
    }
    Z[i] = max_prob_state;
  }
  hiddenStateSeq[totalNumFeatures - 1] = max_prob_state;
  for(i = totalNumFeatures - 1; i > 0; i--){
    int s = Z[i-1];
    hiddenStateSeq[i-1] = T2[s][i];
  }
  
  return CheckMinDuration(features , hiddenStateSeq,           numStates,               allMixtureMeans,            allMixtureVars, 
			  numMixEachState,          totalNumFeatures,        B,                          numElemEachState,
			  T1,                       T2,                      Pi,        mixtureWeight);
}


/******************************************************************************
   Checkforminimumduration() : 
   inputs : Observation Space (pointer to all training feature vectors) same as Sequence of Observations, 
   allMixtureMeans of GMMs - pointer to means of gmms, 
   allMixtureVars of GMMs, (*numStates), numMixEachState, Transition matrix (not using currently), 
   Emission probability matrix Bj(Ot), Initial probabilities (Pi)

   outputs : check for minimum duration and drops any state which has less than MIN_DUR frames 
******************************************************************************/
int* CheckMinDuration(VECTOR_OF_F_VECTORS *features,     int *hiddenStateSeq, int *numStates, 
		      VECTOR_OF_F_VECTORS *allMixtureMeans, VECTOR_OF_F_VECTORS *allMixtureVars, int *numMixEachState, 
		      int totalNumFeatures, float **B, int *numElemEachState, float *T1[], int *T2[], float *Pi,
		      float *mixtureWeight){
  
  FindNumberOfElemInEachState(hiddenStateSeq, numStates, totalNumFeatures,  numElemEachState, Pi);
  // check min duration constraint 
  int              i = 0, j = 0, s = 0, d= 0, mixCount = 0;
  for(s = 0; s < *numStates; s++){
    if( numElemEachState[s] < MIN_DUR ){
      printf("dropping state: %d ....does not contain enough elements...\n", s);
      //drop the state 
      VECTOR_OF_F_VECTORS                  *tempAllMixtureMeans, *tempAllMixtureVars;
      tempAllMixtureMeans   = (VECTOR_OF_F_VECTORS *) calloc(MAX_NUM_MIX * (*numStates) , sizeof(VECTOR_OF_F_VECTORS));
      tempAllMixtureVars    = (VECTOR_OF_F_VECTORS *) calloc(MAX_NUM_MIX * (*numStates) , sizeof(VECTOR_OF_F_VECTORS));
      for(i = 0; i < MAX_NUM_MIX * (*numStates); i++){
	tempAllMixtureMeans[i] = (F_VECTOR *) AllocFVector(DIM);
	tempAllMixtureVars[i]  = (F_VECTOR *) AllocFVector(DIM);
      }
      // COPY ORIGINAL MEANS AND VARS INTO TEMPS
      int totalMix = 0;
      for(i = 0; i < *numStates; i++)
	totalMix += numMixEachState[i];
      
      for(i = 0; i < totalMix; i++){
	for(d = 0; d < DIM; d++){
	  tempAllMixtureMeans[i]->array[d] = allMixtureMeans[i]->array[d];
	  tempAllMixtureMeans[i]->numElements = DIM;
	  tempAllMixtureVars[i]->array[d]  = allMixtureVars[i]->array[d];
	  tempAllMixtureMeans[i]->numElements = DIM;
	}
      }
      
      // copy from temp to original means and vars array
      int mix_s = numMixEachState[s];
      for(j = 0; j < *numStates - 1; j++){
	if(j >= s)
	  numMixEachState[j] = numMixEachState[j+1];
      }
      
      for(j = 0; j < *numStates - 1; j++){
	if(j < s){
	  Pi[j] = (float) log((float)numElemEachState[j]/totalNumFeatures);
	  mixCount = 0;
	  for(i = 0; i < j; i++)
	    mixCount += numMixEachState[i];   
	  for(i = mixCount; i < mixCount + numMixEachState[j]; i++){
	    for(d = 0; d < DIM; d++){
	      allMixtureMeans[i]->array[d] = tempAllMixtureMeans[i]->array[d];
	      allMixtureMeans[i]->numElements = DIM;
	      allMixtureVars[i]->array[d]  = tempAllMixtureVars[i]->array[d];
	      allMixtureMeans[i]->numElements = DIM;
	    }
	  }
	}
	else if(j >= s){
	  mixCount = 0;
	  numElemEachState[j] = numElemEachState[j+1];
	  Pi[j] = (float) log((float)numElemEachState[j]/totalNumFeatures);
	  //printf("Pi[%d] : %f\n", j, Pi[j]);
	  for(i = 0; i < j; i++)
	    mixCount += numMixEachState[i];
	  for(i = mixCount; i < mixCount + numMixEachState[j]; i++){
	    for(d = 0; d < DIM; d++){
	      allMixtureMeans[i]->array[d] = tempAllMixtureMeans[i + mix_s]->array[d];
	      allMixtureMeans[i]->numElements = DIM;
	      allMixtureVars[i]->array[d]  = tempAllMixtureVars[i + mix_s]->array[d];
	      allMixtureMeans[i]->numElements = DIM;
	    }
	  }
	}
      }
      // FREE TEMP VECTORS
      for(i = 0; i < MAX_NUM_MIX * (*numStates); i++){
	free(tempAllMixtureMeans[i]);
	free(tempAllMixtureVars[i]);
      }
      free(tempAllMixtureMeans);
      free(tempAllMixtureVars);
      // change number of states by one
      *numStates = *numStates - 1;
      // first compute new posterior probs for each element
      ComputePosteriorProb(features,    B,      allMixtureMeans,       allMixtureVars,     numStates,       numMixEachState,
			   totalNumFeatures,    mixtureWeight);
      // Now call viterbi alginment again 
      return viterbi(  features, numStates,       allMixtureMeans,      allMixtureVars,      numMixEachState,
		       totalNumFeatures,          B,                    numElemEachState,        T1,           T2,
		       Pi, 		          mixtureWeight);
    }//matches if 
  }
  return hiddenStateSeq;
}
