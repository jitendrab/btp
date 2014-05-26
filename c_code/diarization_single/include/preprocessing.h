#ifndef PREPROCESSING_H
#define PREPROCESSING_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
//#define  ARMA_DONT_USE_WRAPPER
//#include <armadillo>
#include "MyTypes.h"
//#include "viterbi_realign.h"
#include <mlpack/core.hpp>
#include <math.h>
#include <limits>
#include "config_single.h"  //for single gaussian
using namespace::mlpack;
using namespace mlpack::util;
using namespace::mlpack::distribution;
using namespace::arma;

//globals
//mat                                  featuresForClustering; // collect features in this matrix
char                                 fileName[100];
float                                probScaleFactor = 1.0;
int                                  ditherMean = 1;
int                                  varianceNormalize = 1;

int                                  dfs = 0;
int                                  **featureSpace;
int                                  allFeaturesCount = 0;
//function prototypes 
int InitializeGMMs(VECTOR_OF_F_VECTORS *, int , int , int *);

void CalculateDeltaBIC(ESHMM *mdHMM, double **deltaBIC, VECTOR_OF_F_VECTORS *features, int T, int *path);

void ComputePosteriorProb( ESHMM *mdhmm, VECTOR_OF_F_VECTORS *features, int totalFeatures, double **posterior);

void ClusteringAndMerging(VECTOR_OF_F_VECTORS *features, int *path, ESHMM *mdHMM, int totalNumFeatures, double **posterior);

int hmmMergeTwoStates(ESHMM *mdHMM, int *path, int totalNumFeatures, double **deltaBIC, VECTOR_OF_F_VECTORS *feat);

void hmm2MinDurationHMM(ESHMM *mdHMM);

void printHMM(ESHMM *mdHMM);

#endif
