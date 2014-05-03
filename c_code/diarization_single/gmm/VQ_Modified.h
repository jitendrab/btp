#ifndef VQ_H
#define VQ_H
#include "FrontEndDefs.h"
//#include "FrontEndTypes.h"
void InitVQ (VECTOR_OF_F_VECTORS *vfv, 
	     VECTOR_OF_F_VECTORS *clusterMeans, 
	     VECTOR_OF_F_VECTORS *clusterVars,int num_clusters, float seed); 
float ComputeDiscriminant(F_VECTOR *clusterMean, 
			  F_VECTOR *clusterVar, F_VECTOR *fvect, int varianceNormalize);
int DecideWhichCluster(F_VECTOR *fvect, VECTOR_OF_F_VECTORS *clusterMeans,
		       VECTOR_OF_F_VECTORS *clusterVars, int num_clusters,
		       int varianceNormalize);
void ComputeVQ(VECTOR_OF_F_VECTORS *vfv, int num_vectors, VECTOR_OF_F_VECTORS *cluster_means, VECTOR_OF_F_VECTORS *cluster_vars, 
float *cluster_elem_cnt, int num_clusters, int varianceNormalize, 
float ditherMean, int iterations, float seed);
#endif




