/*-------------------------------------------------------------------------
 *  GMM.h - Definition module for GMM.c
 *  Version:	$Name$
 *  Module:	
 *
 *  Purpose:	
 *  See:	
 *
 *  Author:	 (hema@localhost.localdomain)
 *
 *  Created:        Sat 25-May-2002 09:09:52
 *  Last modified:  Wed 07-May-2003 15:51:27 by hema
 *  $Id$
 *
 *  Bugs:	
 *
 *  Change Log:	<Date> <Author>
 *  		<Changes>
 -------------------------------------------------------------------------*/

#ifndef GMM_H
#define GMM_H
void InitGMM (VECTOR_OF_F_VECTORS *vfv, 
	     VECTOR_OF_F_VECTORS *clusterMeans, 
	     VECTOR_OF_F_VECTORS *clusterVars, int numClusters, int seed); 
float ComputeProbability(F_VECTOR *meanVector, F_VECTOR *varVector, 
			 float priorProb, F_VECTOR *fvect, 
			 float probScaleFactor);
int DecideWhichMixture(F_VECTOR *fvect, VECTOR_OF_F_VECTORS *clusterMeans, 
		       VECTOR_OF_F_VECTORS *clusterVars, int numClusters,
		       float *clusterElemCnt, int numVectors,
		       float probScaleFactor);
void ComputeGMM(VECTOR_OF_F_VECTORS *vfv, int numVectors, 
		VECTOR_OF_F_VECTORS *clusterMeans, 
		VECTOR_OF_F_VECTORS *clusterVars, 
		float *clusterElemCnt, int numClusters, 
		int VQIter, int GMMIter, float probScaleFactor, 
                int ditherMean, int varianceNormalize, int seed);
#endif


/*-------------------------------------------------------------------------
 * $Log$
 *
 * Local Variables:
 * time-stamp-active: t
 * time-stamp-line-limit: 20
 * time-stamp-start: "Last modified:[ 	]+"
 * time-stamp-format: "%3a %02d-%3b-%:y %02H:%02M:%02S by %u"
 * time-stamp-end: "$"
 * End:
 *                        End of GMM.h
 -------------------------------------------------------------------------*/
