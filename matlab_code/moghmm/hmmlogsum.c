/* Sum of logarithmic numbers.
 *
 * Implementation based on the paper 'Numerically stable Hidden Markov
 * Model Implementation', T.P. Mann, 2006
 *
 * Dxd, 13-7-2009
 */

#include "mex.h"
#include <math.h>

/* General constants */
double LOG0;


/* Sum two log doubles */ 
double sumlog(double logx, double logy)
{
	double s;
	if (logx==LOG0)
		return logy;
	if (logy==LOG0)
		return logx;
	if (logx>=logy)
		s = logx + log(1.0 + exp(logy-logx));
	else
		s = logy + log(1.0 + exp(logx-logy));
	return s;
}


void mexFunction(int nlhs,mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
	/* Do checks */
	if ((nrhs<1)||(nrhs>2))
		mexErrMsgTxt("Only one or two input arguments: logB = logsum(logA,dim)");
	if (nlhs!=1)
		mexErrMsgTxt("One output argument required: logB = logsum(logA,dim)");
	if (mxGetNumberOfDimensions(prhs[0])!=2)
		mexErrMsgTxt("The first input arguments has to have two dimensions");

	/* Initialize */
	LOG0 = log(0.0);
	long M = mxGetM(prhs[0]);
	long N = mxGetN(prhs[0]);
	long dim;
	double *X = mxGetPr(prhs[0]);
	double *sum;
	int i,j;
	nlhs = 1;

	/* Check if we have to sum the rows */
	if (nrhs<2)
		dim = 1;
	else
		dim = mxGetScalar(prhs[1]);

	/* Now sum vertically or horizontally */
	if (dim==1)
	{
		/* Sum over the rows */
		plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
		sum = mxGetPr(plhs[0]);
		for (j=0; j<N; j++) {
			sum[j] = X[j*M];
			for (i=1; i<M; i++)
				sum[j] = sumlog(sum[j], X[j*M+i]);
		}
	}
	else
	{
		/* Sum over the columns */
		plhs[0] = mxCreateDoubleMatrix(M,1,mxREAL);
		sum = mxGetPr(plhs[0]);
		for (j=0; j<M; j++) {
			sum[j] = X[j];
			for (i=1; i<N; i++)
				sum[j] = sumlog(sum[j], X[j+i*M]);
		}
	}
}
