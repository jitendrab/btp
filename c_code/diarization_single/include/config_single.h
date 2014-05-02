#ifndef CONFIG_SINGLE
#define CONFIG_SINGLE

#define Max_Segs 5000 // max number of voiced speech segments in .scp file 
#define DIM 19
#define MAX_NUM_FEATURES 200000
#define MAX_NUM_STATES 16
#define MAX_NUM_MIX 100
#define INITIAL_NUM_MIX 1 //single gaussian
#define LAMBDA 8
#define INITIAL_STATES 9
#define MIN_DUR 250

typedef struct Gauss{
  F_VECTOR *mean;
  VECTOR_OF_F_VECTORS *cov;
  // functions which will be used for gaussian class in c++
  // void train_gaussian(VECTOR_OF_F_VECTORS *features, F_VECTOR *mean, VECTOR_OF_F_VECTORS *cov);
  // float ComputProbability(F_VECTOR *mean, VECTOR_OF_F_VECTORS *cov, F_VECTOR *fvect, float prior, float probScaleFactor);
  
}gaussian;

typedef struct ESHMM{
  gaussian *states;
  int *numStates;
  VECTOR_OF_F_VECTORS *trans; //transition matrix
  VECTOR_OF_F_VECTORS *prior; // prior probabilities
  int MD;
  // functions which will be used by ESHMM class
  // 3. Build HMM given observation sequence training problem of hmm
  //eshmm* TrainESHMM(eshmm *hmm, VECTOR_OF_F_VECTORS *obserSeq, int *totalNumFeatures, int *stateSeq);
  // 2. compute new optimal state Sequence given the model, observation seq
  // float* OptimalStateSeq(eshmm *hmm, float **posterior, VECTOR_OF_F_VECTORS *features);
  // 1. evaluation problem compute total probability give observation seq
  // float ComuteTotalProb(eshmm *hmm, VECTOR_OF_F_VECTORS *features);
  
}eshmm;

#endif
