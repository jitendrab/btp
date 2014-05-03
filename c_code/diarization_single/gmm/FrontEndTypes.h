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

#ifndef FRONT_END_TYPES
#define FRONT_END_TYPES
#include "stdio.h" 
typedef int INT_TYPE;

typedef float FLOAT_TYPE;

/* typedef struct { */
/*   int numElements; */
/*   float *array; */
  
/* } F_VECTOR; */

/* typedef F_VECTOR* VECTOR_OF_F_VECTORS; */
struct F_VECTOR{
  int numElements;
  float *array; 
 
  /* F_VECTOR(){ */
  /*   numElements = 0; */
  /* } */
  F_VECTOR(int numElems){
    numElements = numElems;
    array = new float[numElems];
  }
  ~F_VECTOR(){
    delete [] array;
  }
};

typedef F_VECTOR* VECTOR_OF_F_VECTORS;

typedef union {
      int ival;
      float fval;
      union uTag *next;
    } u;

typedef struct {
  char  waveFileName[500];
  short *waveform;
  VECTOR_OF_F_VECTORS *melCepstrumCosineTransform,
    *melCepstrumInvTransform;
  VECTOR_OF_F_VECTORS *filterbankWeights;
  short *vU;
  int *dftIndices;
  int fileChanged;
  unsigned int waveType;
  unsigned int windowSize;
  int resGdWindowSize;
  unsigned int fftSize;
  unsigned int fftOrder;
  float preemphasis;
  unsigned int preemphasisDelay;
  unsigned long numSamples;
  unsigned int frameAdvanceSamples;
  unsigned int numCepstrum;
  unsigned int numFilters;
  unsigned int numRegressCoeffts;
  unsigned int numFrames;
  unsigned int numVoicedFrames;
  unsigned int samplingRate;
  unsigned int lpOrder;
  unsigned int zeroOrder;
  unsigned int filterOrder;
  unsigned int minPitch;
  unsigned int maxPitch;
  unsigned int numFormants;
  unsigned int seed;
  unsigned int deltaDifference;
  unsigned int deltaDeltaDifference;
  float filterWarp;
  float trapezoidalRatio;
  float bandwidthScale;
  float gamma;
  float gdPosScale;
  float gdNegScale;
  int gdSmthWinSize;
  unsigned int gdLifterWinSize;
  unsigned int gdRemoveLPhase;
  unsigned int removeMin;
  unsigned int gdSign;
  unsigned int mgdNormalize;
  unsigned int medianOrder;
  unsigned int zeroMean;
  unsigned int varianceNormalize;
  unsigned int featureVarNormalize;
  unsigned int percentFrames;
  unsigned int vad;
  unsigned int perceptualFilterbank;
  unsigned int stGauss;
  unsigned int stGaussWnd;
  unsigned int slopeDCT;
  float varianceFloor;
  float ditherMean;
  float winScaleFactor;
  float probScaleFactor;
  float minFrequency, maxFrequency;
  float thresEnergy;
  float thresZero;
  float thresSpecFlatness;
} ASDF;


typedef struct {
  int numElements;
  int *array;
} I_VECTOR;

typedef I_VECTOR* VECTOR_OF_I_VECTORS;

typedef struct {
  int numColumns;
  int numRows;
  float **array;
} MATRIX;
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
