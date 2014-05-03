/*-------------------------------------------------------------------------
 *  front-end-types.h - A set of definitions used by the speech front-end
 *  Version:	$Name:  $
 *  Module:	
 *
 *  Purpose:	
 *  See:	
 *
 *  Author:	Hema A Murthy (hema@bhairavi.iitm.ernet.in)
 *
 *  Created:        Mon 10-Sep-2001 15:33:45
 *  Last modified:  Fri 25-Mar-2011 14:02:42 by hema
 *  $Id: front-end-types.h,v 1.1 2001/10/25 12:28:05 hema Exp $
 *
 *  Bugs:	
 *
 *  Change Log:	<Date> <Author>
 *  		<Changes>
 -------------------------------------------------------------------------*/

#ifndef MY_TYPES
#define MY_TYPES
#include <stdio.h>
using namespace::std;
struct F_VECTOR{
 public:
  int numElements;
  float *array;
  
  F_VECTOR();
  F_VECTOR(int numElems){
    numElements = numElems;
    array = new float[numElems];
  }
  ~F_VECTOR(){
    delete [] array;
  }
};

typedef F_VECTOR* VECTOR_OF_F_VECTORS;

class Gaussian{
 public:
  F_VECTOR *mean;
  VECTOR_OF_F_VECTORS *cov;
  //constructor
 public:
  Gaussian(int DIM){    
    mean = new F_VECTOR(DIM);
    VECTOR_OF_F_VECTORS *cov = new F_VECTOR *[DIM];
    int i = 0;
    for(i = 0; i < DIM; i++)
      cov[i] = new F_VECTOR(DIM);    
  }
  //destructor
  ~Gaussian(){
    int numElems = mean->numElements;
    int i = 0;
    delete mean;
    for(i = 0; i < numElems; i++)
      delete cov[i];
    delete [] cov;
  }
  
  // Functions which will be used for gaussian class in c++
  // void train_gaussian(VECTOR_OF_F_VECTORS *features, F_VECTOR *mean, VECTOR_OF_F_VECTORS *cov);
  // float ComputProbability(F_VECTOR *mean, VECTOR_OF_F_VECTORS *cov, F_VECTOR *fvect, float prior, float probScaleFactor);  
};

class ESHMM{
 public:
  Gaussian **HMMstates;
  int numStates;
  VECTOR_OF_F_VECTORS *trans; //transition matrix
  VECTOR_OF_F_VECTORS *prior; // prior probabilities
  int MD;
  int rowsTrans, rowsPrior;
 public:
  ESHMM();
  ESHMM(int numStates, int DIM, int rowsTrans, int colsTrans, int rowsPrior, int colsPrior){
    numStates = numStates;
    rowsTrans = rowsTrans;
    rowsPrior = rowsPrior;
    // HMMstates is an array of pointers pointing to numStates Gaussian objects
    // first allocate numStates pointers to Gaussian object and then allocate Gaussian object for each pointer
    HMMstates = new Gaussian *[numStates];
    int i = 0;
    for(i = 0; i < numStates; i++)
      HMMstates[0] = new Gaussian(DIM);
    
    trans = new F_VECTOR *[rowsTrans];
    prior = new F_VECTOR *[rowsPrior];    
    for(i = 0; i < rowsTrans; i++)
      trans[i] = new F_VECTOR(colsTrans);    
    for(i = 0; i < rowsPrior; i++)
      prior[i] = new F_VECTOR(colsPrior);
  }
    
  ~ESHMM(){
    int i = 0;
    for(i = 0; i < rowsTrans; i++)
      delete [] trans[i];
    for(i = 0; i < rowsPrior; i++)
      delete [] prior[i];
    for(i = 0; i < numStates; i++)
      delete [] HMMstates[i];
    delete [] prior;
    delete [] trans;
    delete [] HMMstates;
  }
  // functions which will be used by ESHMM class
  // 3. Build HMM given observation sequence training problem of hmm
  //eshmm* TrainESHMM(eshmm *hmm, VECTOR_OF_F_VECTORS *obserSeq, int *totalNumFeatures, int *stateSeq);
  // 2. compute new optimal state Sequence given the model, observation seq
  // float* OptimalStateSeq(eshmm *hmm, float **posterior, VECTOR_OF_F_VECTORS *features);
  // 1. evaluation problem compute total probability give observation seq
  // float ComuteTotalProb(eshmm *hmm, VECTOR_OF_F_VECTORS *features);
  
};

/* typedef struct matrix{ */
/*   int numColumns; */
/*   int numRows; */
/*   float **array; */
/* } MATRIX; */

#endif


/*-------------------------------------------------------------------------
 * $Log: front-end-types.h,v $
 * Revision 1.1  2001/10/25 12:28:05  hema
 * Initial revision
 *
 *
 * Local Variables:
 * time-stamp-active: t
 * time-stamp-line-limit: 20
 * time-stamp-start: "Last modified:[ 	]+"
 * time-stamp-format: "%3a %02d-%3b-%:y %02H:%02M:%02S by %u"
 * time-stamp-end: "$"
 * End:
 *                        End of front-end-types.h
 -------------------------------------------------------------------------*/
