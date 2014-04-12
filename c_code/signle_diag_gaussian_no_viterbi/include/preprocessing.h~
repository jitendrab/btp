#ifndef PREPROCESSING_H
#define PREPROCESSING_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FrontEndDefs.h"
#include "FrontEndTypes.h"
#include "VQ_Modified.h"
#include "InitAsdf.h"
#include "DspLibrary.h"
#include "GMM.h"
#include "VQ_Modified.h"
//#include "viterbi_realign.h"
#include <math.h>

#include "config_single.h"  //for single gaussian
//globals
VECTOR_OF_F_VECTORS                  *mixtureMeans, *mixtureVars;
F_VECTOR                             *mean;
F_VECTOR                             *var;
F_VECTOR                             *x;
VECTOR_OF_F_VECTORS                  *featuresForClustering;
char                                 fileName[100];
float                                probScaleFactor = 1.0;
int                                  ditherMean = 1;
int                                  varianceNormalize = 1;
int                                  featureSpace[MAX_NUM_FEATURES][3];
int                                  allFeaturesCount = 0;
//function prototypes 

int InitializeGMMs(VECTOR_OF_F_VECTORS *, int , int , int *);

void ComputePosteriorProb( VECTOR_OF_F_VECTORS *features, float **posterior, VECTOR_OF_F_VECTORS *allMixtureMeans, 
			   VECTOR_OF_F_VECTORS *allMixtureVars, int *numStates,  int );

int ClusteringAndMerging(VECTOR_OF_F_VECTORS *features,          VECTOR_OF_F_VECTORS *allMixtureMeans,    
			  VECTOR_OF_F_VECTORS *allMixtureVars,               int *numStates,   
			  int totalNumFeatures,         float **posterior,   int *numElemEachState,
			  float *Pi);

//void CalculateBIC(float **BIC, VECTOR_OF_F_VECTORS *,    int *numStates,  float **posterior, int totalNumFeatures, int *stateSeq, int *numMixEachState,               int *numElemEachState);

void FindNumberOfElemInEachState(int *newStateSeq, int *numStates, int totalNumFeatures,  int *numElemEachState, float *Pi);

/*void MergeTwoStates( int s_i,  int s_j  , VECTOR_OF_F_VECTORS *allMixtureMeans,     
		     VECTOR_OF_F_VECTORS *allMixtureVars,      int *numStates,    int *numMixEachState,  
		     int *numElemEachState,                    float *Pi,         int );*/

void CalculateBIC( float **deltaBIC,   VECTOR_OF_F_VECTORS *features,    VECTOR_OF_F_VECTORS *allMixtureMeans, 
		   VECTOR_OF_F_VECTORS *allMixtureVars,                  int *numStates,   int totalNumFeature, int *stateSeq,  
		   int *numMixEachState,               int *numElemEachState, float **posterior);

void BIC_Modified( float **deltaBIC,   VECTOR_OF_F_VECTORS *features,    VECTOR_OF_F_VECTORS *allMixtureMeans, 
		   VECTOR_OF_F_VECTORS *allMixtureVars,                  int *numStates,   int totalNumFeature, int *stateSeq,  
		   int *numElemEachState);

void MergeTwoStates( int s_i,  int s_j  , VECTOR_OF_F_VECTORS *allMixtureMeans,     
		     VECTOR_OF_F_VECTORS *allMixtureVars,      int *numStates,    int *numMixEachState,  
			     int *numElemEachState,                    float *Pi,          int totalNumFeatures,
		     VECTOR_OF_F_VECTORS *features, int *stateSeq, float *mixtureWeight);

void MergeTwoStatesModified( int s_i,  int s_j  , VECTOR_OF_F_VECTORS *allMixtureMeans,     
		     VECTOR_OF_F_VECTORS *allMixtureVars,      int *numStates,   
			     int *numElemEachState,                    float *Pi,          int totalNumFeatures,
			     VECTOR_OF_F_VECTORS *features, int *stateSeq);

void writeRTTMFile(int *stateSeq, int *numStates, int totalNumFeatures, int *numElemEachState);

void writePlotFile(int totalNumFeatures, int *numStates, int *);

void PrintAllDetails(int *numStates, int *numMixEachState, int *numElemEachState, 
		     float **mixtureElemCount, VECTOR_OF_F_VECTORS *allMixtureMeans, float *mixtureWeight);

#endif
