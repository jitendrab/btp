/*-------------------------------------------------------------------------
 *  DspLibrary.c - A Library of Signal Processing Programs
 *  Version:	$Name:  $
 *  Module:	
 *
 *  Purpose:	
 *  See:	
 *
 *  Author:	Hema A Murthy (hema@bhairavi.iitm.ernet.in)
 *
 *  Created:        Some Time in 1996
 *  Last modified:  Thu 11-Mar-2004 15:52:25 by hema
 *  $Id: DspLibrary.c,v 1.2 2004/12/23 08:28:52 hema Exp hema $
 *
 *  Bugs:	
 *
 *  Change Log:	<Date> <Author>
 *  		<Changes>
 -------------------------------------------------------------------------*/
#include "stdio.h"
#include "math.h"
#include "malloc.h"
#include "fe/constants.h"
#include "fe/FrontEndDefs.h"
#include "fe/FrontEndTypes.h"
#include "fe/DspLibrary.h"

float HanW(int i,int npts) { 
 return(0.5+0.5*cos(PI*((i)-1)/((npts)-1)));
}

float HanDw(int i,int npts){ 
  return((float)(0.5-0.5*cos((double) PI2*((i)-1)/((npts)-1))));
}

float HamW(int i,int npts) {
  return(0.54+0.46*cos((double) PI*((i)-1)/((npts)-1)));
}

float HamDw(int i,int npts) {
  return((float) (0.54-0.46*cos((double) PI2*((i)-1)/((npts)-1))) );
}

float BartW(int i,int npts) {
  return (((float)((double) (npts)-(i)))/((npts)-1));
}

/*--------------------------------------------------------------------------- 
  Window applies a specified window to a signal[]  

  inputs - 
            signal - array of floating point no's
            npts  -  size of signal array
            hw - 'M' - Hamming window
                 'N' - Hann window
                 'B' - Bartlett window
                 'G' - Gaussian window 
                 'R' - Rectangular window 
     All other values invalid and taken as 'rectangular' by default.

            sw - 'S' - single sided
                 'D' - double sided
                       all other values invalid   
            gausmin - for Gaussian window only: the least value at the 
                  window ends i.e 1st and last values.The parameter 
                  'min' is ignored for other windows.
 --------------------------------------------------------------------------*/

 void Window(float *sig,int npts,char hw,char sw,float gausmin)
 {
   int i;
   float midpoint,sigma2,tempval,temp;

   if ((hw=='M') && (sw=='S')) /* Single sided ham Window */
      for(i=1; i<=npts; i++ ) sig[i]=sig[i]*HamW(i,npts);
   
   else if ((hw=='M') && (sw=='D')) /* Double sided ham Window */
      for(i=1; i<=npts; i++ ) sig[i]=sig[i]*HamDw(i,npts); 

   else if ((hw=='N') && (sw=='S')) /* Single sided han Window */
      for(i=1; i<=npts; i++ ) sig[i]=sig[i]*HanW(i,npts);

   else if ((hw=='N') && (sw=='D')) /* Double sided han Window */
      for(i=1; i<=npts; i++ ) sig[i]=sig[i]*HanDw(i,npts);
   
   else if ((hw=='B') && (sw=='S')) /* Single sided Bartlett Window */
      for(i=1;i<=npts;i++) sig[i]=sig[i]*BartW(i,npts);

   else if ((hw=='B') && (sw=='D')){ /* Double sided Bartlett Window */
      midpoint=((float)(npts)+1.0)/2.0;
      for(i=1;i<=npts;i++){
        tempval=(2.0*(i-1))/(float)(npts-1);
        sig[i]*=(((float)(i))<midpoint?tempval:(2.0-tempval));
      }
   } 
   else if (hw == 'G'){
     if(gausmin < 1.0){
        tempval = sqrt((double)(-2.0*log(gausmin)));
        if(sw =='S'){
          sigma2 = ((float)npts-1.0)/tempval; /* Standard Dev. */
          sigma2 = -2.0*sigma2*sigma2; 
          for(i = 1; i <= npts; i++){
             temp = exp((double)((i-1)*(i-1)/sigma2));
             sig[i]*= temp;
          }
        }
        else if(sw == 'D'){
          midpoint = ((float)(npts)+1.0)/2.0;
          sigma2 = (midpoint-1.0)/tempval;  /* Standard Dev. */
          sigma2 = -2.0*sigma2*sigma2;
          for(i = 1;i <= npts; i++){
            temp = exp((double)(((float)i-midpoint)*((float)i-midpoint)/sigma2));
            sig[i]*= temp;
          }
        }
        else printf("\n ERROR from Window: Bad Win. Symmetry.\n");
     }
     else printf("\n Error from Gaus.win.:'min' cannot be < 1 .");
   }
   else if (hw!='R') {
      printf("Warning from subroutine Window:\n"); 
      printf("\t\t error in Window specification\n");
      printf("\t\t Rectangular Window is used.\n");
   }
 }


/*--------------------------------------------------------------------------
    The subroutine TrapeziumWindow() Windows the passed signal with
    a Window which is rectangular in the middle and has a taper
    at the ends of the specified type,ex:Hamming,Hann,Bartlett etc.
    Inputs -
           sig[] - (float)array of signal to be Windowed.
           npts - (int) total Window length = winSep + 2*winLen
           taperWin - (char) type of Window used as taper at ends.
           winSep - (int)the length of the central portion in which it
                     is rectangular.
    Output -
           sig[] - the passed signal array is Windowed and passed
                   back.
    Note: All arrays are indexed from 1 and not zero.
 --------------------------------------------------------------------------*/

 void TrapeziumWindow(float sig[],int npts,char taperWin,int winSep)
 {
    int n,dummy,endWinLen;
    float temp[256+1];

    /* get taper Window(one sided) into temp[]. */
    endWinLen=(npts-winSep)/2;  /* this is length of endWindows. */
    for(n=1;n<=endWinLen;temp[n++]=1.0);
    Window(temp,endWinLen,taperWin,'S',0.01);

    for(n=1;n<=endWinLen;n++)
       sig[n]*=temp[endWinLen-n+1];

    dummy=endWinLen+winSep;
    for(n=dummy+1;n<=npts;n++)
       sig[n]*=temp[n-dummy];
 }

/*--------------------------------------------------------------------------
  Imin returns the location of the minimum of an array signal of 
         size npts.                          
 -------------------------------------------------------------------------*/

 int Imin(float sig[],int npts)
 {
   int i,lmin;
   lmin = 1;
   for ( i=2; i<=npts; i++ ) 
     if (fabs(sig[i]) < fabs(sig[lmin])) lmin = i;
   return (lmin);
 }

/*--------------------------------------------------------------------------
    Imax returns the location of the maximum of 
    an array signal of size npts. 
 -------------------------------------------------------------------------*/
 int Imax(float sig[],int npts)
 {
   int i,lmax;
   lmax = 1;
   for ( i=2; i<=npts; i++ ) 
     if (fabs(sig[i]) > fabs(sig[lmax])) lmax = i;
   return (lmax);
 }
/*--------------------------------------------------------------------------
  Imin returns the location of the minimum of an array signal of 
         size npts(indices start from 0)                          
 -------------------------------------------------------------------------*/

 int Imin0(float sig[],int npts)
 {
   int i,lmin;
   lmin = 0;
   for ( i=1; i<npts; i++ ) 
     if (fabs(sig[i]) < fabs(sig[lmin])) lmin = i;
   return (lmin);
 }

/*--------------------------------------------------------------------------
    Imax returns the location of the maximum of 
    an array signal of size npts (indices start from 0) 
 -------------------------------------------------------------------------*/
 int Imax0(float sig[],int npts)
 {
   int i,lmax;
   lmax = 0;
   for ( i=1; i<npts; i++ ) 
     if (fabs(sig[i]) > fabs(sig[lmax])) lmax = i;
   return (lmax);
 }
/*--------------------------------------------------------------------------
  Imin returns the location of the minimum of an array signal of 
         size npts.                          
 -------------------------------------------------------------------------*/

 int IminActual(float sig[],int npts)
 {
   int i,lmin;
   lmin = 1;
   for ( i=2; i<=npts; i++ ) 
     if (sig[i] < sig[lmin]) lmin = i;
   return (lmin);
 }



/*--------------------------------------------------------------------------
    Imax returns the location of the maximum of 
    an array signal of size npts. 
 -------------------------------------------------------------------------*/
 int ImaxActual(float sig[],int npts)
 {
   int i,lmax;
   lmax = 1;
   for ( i=2; i<=npts; i++ ) 
     if (sig[i] > sig[lmax]) lmax = i;
   return (lmax);
 }
/*--------------------------------------------------------------------------
  Imin returns the location of the minimum of an array signal of 
         size npts.                          
 -------------------------------------------------------------------------*/

 int Imin0Actual(float sig[],int npts)
 {
   int i,lmin;
   lmin = 0;
   for ( i=1; i<npts; i++ ) 
     if (sig[i] < sig[lmin]) lmin = i;
   return (lmin);
 }

/*--------------------------------------------------------------------------
    FindIndex checks whether a given element is found in a 
    an array of ints of size npts
 -------------------------------------------------------------------------*/
 int FindIndex(int *array,int npts, int index)
 {
   int flag = 0;
   int i = 0;
   while ((i <= npts) && (!flag)) 
     if (array[i] == index)
       flag = 1;
     else 
       i++;
   return(flag);
 }
/*--------------------------------------------------------------------------
    FindMatch checks whether a given element already exists 
    an array of ints of size npts
 -------------------------------------------------------------------------*/
 int FindMatch(VECTOR_OF_F_VECTORS *vfv, int numVectors, int *array,int npts, int index)
 {
   F_VECTOR *fvect1, *fvect2;
   fvect1 = vfv[index];
   int flag = 0;
   int i = 0, k;
   float sum;
   while ((i <= npts) && (!flag)){
     fvect2 = vfv[i];
     k = 0;
     sum = 0;
     while (k < fvect1->numElements) {
       sum = sum + fvect1->array[k] - fvect2->array[k];
       k++;
     }
     if (sum == 0) flag = 1;
   i++;
   }
   return(flag);
 }

/*--------------------------------------------------------------------------
    Imax returns the location of the maximum of 
    an array signal of size npts. 
 -------------------------------------------------------------------------*/
 int Imax0Actual(float sig[],int npts)
 {
   int i,lmax;
   lmax = 0;
   for ( i=1; i<npts; i++ ) 
     if (sig[i] > sig[lmax]) lmax = i;
   return (lmax);
 }

/*-------------------------------------------------------------------------
median returns the location of the middle value 
in an array signal of size npts
--------------------------------------------------------------------------*/
 float Median(float sig[],int npts)
 {
   static           float *sigcopy;
  float             temp;
  int               min,i,j;
  static int        flag = 0;

  if (flag == 0) {
    sigcopy = (float *) AllocFloatArray(sigcopy, npts+1);
    flag = 1;
  }
  for (i = 1; i <= npts; i++)
    sigcopy[i] = sig[i];
  for (i = 1; i <= npts; i++) {
    min = i;
    for (j = i+1; j <= npts; j++)
      if (sigcopy[j] < sigcopy[min]) 
        min = j;
    temp = sigcopy[i];
    sigcopy[i] = sigcopy[min];
    sigcopy[min] = temp;
  }
 return(sigcopy[(npts+1)/2]);
 }

/*------------------------------------------------------------------------- 
   FindPeaks finds the peaks in the given array of data.
   Inputs-
           Spectrum[] - array of length nfft
           nfft - no of points in the input array

   Output-
            randPhs[]  -  array of nfft integers

 Note: This routine can handle a data array of 1024 samples only, as
       the 'slopes[]' array in the routine has only this length.
 ------------------------------------------------------------------------*/
 void FindPeaks(float *Spectrum,int *randPhs,int npts)
 {
   int  i,k;
   randPhs[1] = 0;
   randPhs[npts] = 0; 
   for (i = 2; i < npts; i++) {
      if ((Spectrum[i-1] < Spectrum[i] )&&(Spectrum[i] > Spectrum[i+1]))  {
       randPhs[i] = 1;
      }
      else randPhs[i] = 0; 
      /*  printf("randPhs %d = %d", i, randPhs[i]);
	  fflush(stdout); */
   }
 }

/*-------------------------------------------------------------------------
 *  AllocFloatArray -- Allocates an array of floats
 *    Args:	Array, size of array
 *    Returns:	allocated array
 *    Bugs:	
 * -------------------------------------------------------------------------*/

float * AllocFloatArray(float *array, int npts)
{
  array = (float *) calloc (npts, sizeof(float));
  if (array == NULL) {
    printf("unable to allocate Float array \n");
    exit(-1);
  }
  return(array);
}	/*  End of AllocFloatArray */	

/*-------------------------------------------------------------------------
 *  AllocIntArray -- Allocates an array of Ints
 *    Args:	Array, size of array
 *    Returns:	allocated array
 *    Bugs:	
 * -------------------------------------------------------------------------*/

int *AllocIntArray(int *array, int npts)
{
  array = (int *) calloc(npts,sizeof(int));
  if (array == NULL) {
    printf("unable to allocate Int array \n");
    exit(-1);
  }
  return(array);
}	/*  End of AllocIntArray  */


/* Global definitions for fft computation  */
static     int *iBit;
static     float *twiddleReal, *twiddleImag;

/*--------------------------------------------------------------------------
  Cstore computes the twiddle factors used in FFT computation.
          The twiddle factors are stored in the global arrays
          IBIT, twiddleReal, twiddleImag.                                    
 -------------------------------------------------------------------------*/

 void Cstore(int n)
 /*  int n;*/        /* FFT order */
 {
   int  nv2,nm1,ix,ix1,j,i,k;
   float pi2byn;

   iBit = (int *) AllocIntArray(iBit,n+1);
   twiddleReal = (float *) AllocFloatArray(twiddleReal,n/2+1);
   twiddleImag = (float *) AllocFloatArray(twiddleImag,n/2+1);
   nv2 = n/2;
   nm1 = n-1;
   iBit[1] = 1;
   iBit[n] = n;
   ix = 0;
   for (i=2; i <= nm1; i++){
     j = 0;
     k = nv2;
     ix1 = ix;     
     while (ix1 >= k) { j = j+k; ix1 = ix1-k; k = k/2; };
     ix = ix + k - j;
     iBit[i] = ix + 1;
   };
   pi2byn = (float)(8.0*atan((double)1.0)/(double)n);
   for (i=1; i <= nv2; i++) {
     k = i-1;
     twiddleReal[i] = (float)cos((double)(pi2byn * k));
     twiddleImag[i] = (float)sin((double)(pi2byn * k));
   }
 }


/* ----------------------------------------------------------------------------


Cfft computes the FT of a complex signal.
        inputs - 
                a - complex signal of length n
                  n - FFT order
                m - m such that 2**m = n
                nsign -  -1  forward
                          1  inverse
        
        outputs - 
                 b - complex array of length n            

-----------------------------------------------------------------------------*/

void FFTReal(a,b,m,n,nsign)
     int                   m,n,nsign;
     complex               a[],b[];
{
  int                       nv2,nm1,i,j,ip,k,le,le1,le2,l;
  static                    float twoPowerM;
  static                    float log2; 
  static                    int flag = 0;
  complex                   u,t;

  if ((int)pow(2,m)!=n){
    printf("ERROR from Cfft: 2**m != n\n");
    exit(1);
  }
  if (flag == 0) {
    log2 = log((double)2.0);
    flag = 1;
  } 
  nv2 = n/2;
  nm1 = n-1;
  for ( i=1; i<=n; i++ ) b[iBit[i]] = a[i]; 
  
  for ( i=1; i<=n/2; i+=2 )
    {
      ip = i+1;
      t = b[ip];
      csub( b[ip],b[i],t ); /* b[ip] = b[i] - t  */
      b[iBit[ip]].im = -b[ip].im;
      cadd( b[i],b[i],t );  /* b[i] = b[i] + t   */
      b[iBit[i]].im = -b[i].im;
    };
  
  for( i=1; i<=n/2; i+=4 )
    {
      ip = i+2;
      t = b[ip];
      csub( b[ip],b[i],t ); /* b[ip] = b[i] - t  */
      b[iBit[ip]].im = -b[ip].im;
      cadd( b[i],b[i],t );  /* b[i] = b[i] + t   */
      b[iBit[i]].im = -b[i].im;
   };

 for( i=2; i<=n/2; i+=4 )
    {
     ip = i+2;
     t.re = -nsign * b[ip].im;
     t.im =  nsign * b[ip].re;
     csub( b[ip],b[i],t ); /* b[ip] = b[i] - t  */
      b[iBit[ip]].im = -b[ip].im;
     cadd( b[i],b[i],t );  /* b[i] = b[i] + t   */
      b[iBit[i]].im = -b[i].im;
    };

 for( l=3; l<=m; l++ )
    {
     le2 = (int) (exp(log((double)2.0)*(m-l)) + 
		  (double)0.5);  /* le2 = 2**(m-l) */
     le = (int) (exp(log2*(double)l)+(double)0.5);   /* le = 2**l */
     le1 = le/2;
     for ( j=1; j<=le1; j++ )
        {
         k = (j-1)*le2+1;
         u.re = twiddleReal[k];
         u.im = nsign*twiddleImag[k];
         for ( i=j; i<=n/2; i+=le )
            {
             ip = i+le1;
             cmul(t,b[ip],u);   /*  t = b[ip]*u  */
             csub( b[ip],b[i],t ); /* b[ip] = b[i] - t  */
	     b[iBit[ip]].im = -b[ip].im;
             cadd( b[i],b[i],t );  /* b[i] = b[i] + t   */
             b[iBit[i]].im = -b[i].im;

            };
        };
    };
 if(nsign==1) for ( i=1; i<=n; i++ ) { b[i].re=b[i].re/(float)n;
                                        b[i].im=b[i].im/(float)n; };
}


/* ----------------------------------------------------------------------------


Cfft computes the FT of a complex signal.
        inputs - 
                a - complex signal of length n
                  n - FFT order
                m - m such that 2**m = n
                nsign -  -1  forward
                          1  inverse
        
        outputs - 
                 b - complex array of length n            

-----------------------------------------------------------------------------*/

void Cfft(a,b,m,n,nsign)
     int                   m,n,nsign;
     complex               a[],b[];
{
  int                       nv2,nm1,i,j,ip,k,le,le1,le2,l;
  static                    float log2; 
  static                    int flag = 0;
  complex                   u,t;

  if ((int)pow(2,m)!=n){
    printf("ERROR from Cfft: 2**m != n\n");
    exit(1);
  }
  if (flag == 0) {
    log2 = log((double)2.0);
    flag = 1;
  } 
  nv2 = n/2;
  nm1 = n-1;
  for ( i=1; i<=n; i++ ) b[iBit[i]] = a[i]; 
  
  for ( i=1; i<=n; i+=2 )
    {
      ip = i+1;
      t = b[ip];
      csub( b[ip],b[i],t ); /* b[ip] = b[i] - t  */
      cadd( b[i],b[i],t );  /* b[i] = b[i] + t   */
    };
  
  for( i=1; i<=n; i+=4 )
    {
      ip = i+2;
      t = b[ip];
      csub( b[ip],b[i],t ); /* b[ip] = b[i] - t  */
      cadd( b[i],b[i],t );  /* b[i] = b[i] + t   */
    };

 for( i=2; i<=n; i+=4 )
    {
     ip = i+2;
     t.re = -nsign * b[ip].im;
     t.im =  nsign * b[ip].re;
     csub( b[ip],b[i],t ); /* b[ip] = b[i] - t  */
     cadd( b[i],b[i],t );  /* b[i] = b[i] + t   */
    };

 for( l=3; l<=m; l++ )
    {
     le2 = (int) (exp(log((double)2.0)*(m-l)) + 
		  (double)0.5);  /* le2 = 2**(m-l) */
     le = (int) (exp(log2*(double)l)+(double)0.5);   /* le = 2**l */
     le1 = le/2;
     for ( j=1; j<=le1; j++ )
        {
         k = (j-1)*le2+1;
         u.re = twiddleReal[k];
         u.im = nsign*twiddleImag[k];
         for ( i=j; i<=n; i+=le )
            {
             ip = i+le1;
             cmul(t,b[ip],u);   /*  t = b[ip]*u  */
             csub( b[ip],b[i],t ); /* b[ip] = b[i] - t  */
             cadd( b[i],b[i],t );  /* b[i] = b[i] + t   */
            };
        };
    };
 if(nsign==1) for ( i=1; i<=n; i++ ) { b[i].re=b[i].re/(float)n;
                                        b[i].im=b[i].im/(float)n; };
}



/*---------------------------------------------------------------------------  

    Rfft computes the FT of a real signal   
    inputs - 
             signal - real array of length nfft
             mfft  -    2**mfft = nfft
             nfft  -  order of fft
             nsign -  -1  forward
		       1  inverse
    outputs -
	     ax - array of real part of  FT
	     ay - array of imaginary part of FT             

-----------------------------------------------------------------------------*/


void Rfft(sig,ax,ay,mfft,nfft,nsign)
 float  sig[], ax[], ay[];
 int  mfft,nfft,nsign;
{
  static                complex *a, *b;
  int i;
  static                int flag = 0;
 if (flag == 0) {
   a = (complex *) malloc((nfft+1)*sizeof(complex));
   b = (complex *) malloc((nfft+1)*sizeof(complex));
   if ((a == NULL) || (b == NULL)) {
     printf("unable to allocate complex array \n");
     exit(-1);
   }
   flag = 1;
 }
 for (i=1; i<=nfft; i++ ) { a[i].re = sig[i];
                            a[i].im = 0.0;       };
 
 FFTReal(a,b,mfft,nfft,nsign);
 //Cfft(a,b,mfft,nfft,nsign);
 for ( i=1; i<=nfft; i++ ) { ax[i] = b[i].re; ay[i] = b[i].im; };
}


/*---------------------------------------------------------------------------

  
    SpectrumReal  computes the magnitude and phase from real and imaginary
           parts of FT
    inputs - 
              ax - array of length nfft - real part
              ay - array of length nfft - imaginary part
    outputs -
              amag - array of length nfft  -  magnitude Spectrum
             phase - "        "            -  phase Spectrum

----------------------------------------------------------------------------*/

void SpectrumReal(nfft,ax,ay,amag,phase)
 int             nfft;
 float           ax[],ay[],amag[],phase[];
{
 int             i;
 for ( i = 1; i <= nfft; i++ )
    {
     amag[i] = sqrt(ax[i]*ax[i] + ay[i]*ay[i]);
     if ((ay[i]==ax[i])&&(ax[i]==0.0))
       phase[i] = atan2((double)ay[i],1.0);
     else
       phase[i] = atan2((double)ay[i],(double)ax[i]);
    }
}


/*---------------------------------------------------------------------------  

  SpectrumComplex computes the magnitude and phase of a complex FT   

----------------------------------------------------------------------------*/

void SpectrumComplex(nfft,csig,amag,phase)
 int                nfft;
 complex            csig[];
 float              amag[],phase[];
{
 int                i;
 for (i=1; i<=nfft; i++ ) 
    {
     amag[i] = csig[i].re*csig[i].re + csig[i].im*csig[i].im;
     if (csig[i].re == 0.0) phase[i] = PI/2.0;
     else phase[i] = atan2(csig[i].im,csig[i].re);
    }
}
/*--------------------------------------------------------------------------
  FrameCompCepstrum : Computes the cepstrum of a given frame of speech

  inputs :
     signal, numPts : signal is an array of numPts
     numCepstrum    : number of cepstral coeffts reqd.
     mfft, nfft      : 2**mfft = nfft
     
  outputs :
     cepstrum        : an array of numCepstrum cepstral coeffts

---------------------------------------------------------------------------*/
float *FrameCompCepstrum(float *signal, int numPts, float *cepstrum, int numCepstrum, int mfft, int nfft)
{
 float resAx[1025],resAy[1025],amag[1025],phase[1025];
 int i;
 for (i = numPts+1; i <= nfft; i++)
   signal[i] = 0;
 Rfft(signal, resAx, resAy, mfft, nfft, -1);
 SpectrumReal(nfft, resAx, resAy, amag, phase);
 for ( i=1; i<=nfft; i++ ) amag[i] = log(amag[i]);
 Rfft(amag,resAx,resAy,mfft,nfft,1);
  for ( i=1; i<=numCepstrum; i++ )
     cepstrum[i-1] = resAx[i+1];
  return(cepstrum);
}

/******************************************************************************
   AllocFVector : Allocate space for a variable of type F_VECTOR with an 
   array size of numPts.
   inputs : npts
   outputs : an F_VECTOR of size npts
******************************************************************************/

F_VECTOR *AllocFVector(int npts) {
   F_VECTOR *fVect;
   fVect = (F_VECTOR *) malloc(1*sizeof(F_VECTOR));
   if (fVect == NULL) {
     printf("unable to allocate FVector \n");
     exit(-1);
   }
   fVect->numElements = npts;
   fVect->array = (float *) calloc(npts, sizeof(float));
   if (fVect->array == NULL) {
     printf("unable to allocate FVector array \n");
     exit(-1);
   } else
    return(fVect);
 }
/*------------------------------------------------------------------------
  LinearVectorDifference : Computers the vectorial difference between two F_VECTORs
  inputs :
    fvect1,fvect2 : two F_VECTORs of identical size
  outputs :
    fvect         : vectorial difference fvect1 - fvect2
-------------------------------------------------------------------------*/   
void LinearVectorDifference (F_VECTOR *fvect1, F_VECTOR *fvect2, F_VECTOR *fvect) {
  int                numElements;
  int                i;

  numElements = fvect1->numElements;
  for (i = 0; i < numElements; i++)
    fvect->array[i] = fvect1->array[i] - fvect2->array[i];
}
/*------------------------------------------------------------------------
  LinearVectorAddition : Computers the vectorial addition between two F_VECTORs
  inputs :
    fvect1,fvect2 : two F_VECTORs of identical size
  outputs :
    fvect         : vectorial addition fvect1 + fvect2
-------------------------------------------------------------------------*/   
void LinearVectorAddition (F_VECTOR *fvect1, F_VECTOR *fvect2, F_VECTOR *fvect) {
  int                numElements;
  int                i;

  numElements = fvect1->numElements;
  for (i = 0; i < numElements; i++)
    fvect->array[i] = fvect1->array[i] + fvect2->array[i];
}
/*------------------------------------------------------------------------
  LinearVectorScalarDivide : Divides a vector by a scalar constant
  inputs :
    scalar        : float
    fvect1        : F_VECTOR 
  outputs :
    fvect         : result fvect1/constant
-------------------------------------------------------------------------*/   
void LinearVectorScalarDivide (float scalar, F_VECTOR *fvect1, F_VECTOR *fvect) {
  int                numElements;
  int                i;

  numElements = fvect1->numElements;
  for (i = 0; i < numElements; i++)
    fvect->array[i] = fvect1->array[i]/scalar;
}

/*******************************************************************************
*	Durbin computes the Auto-Regressive coefficients coef of order m and the
*	reflection coefficients refCoef also of order m 
*       Inputs -
*	autoCorr - Autocorrelation coefficients
        m        - order of the prediction filter
        Outputs -
*	resEnergy - residual energy used by LpAnal for gain computation
*       coef      - AR coefficients
*       refCoef   - Reflection Coefficients
******************************************************************************/
void Durbin(double *autoCorr,int m,float *coef,float *refCoef,float *resEnergy){
  static float 		*b;
  double 		sum;
  int     	        i,i1,j,imj,ij1;
  double                Ecopy;
  static int            flag = 0;
  if (flag == 0) { 
    b = (float *) AllocFloatArray(b,(m+2));
    flag = 1;
  }
  Ecopy = autoCorr[1];
  refCoef[1] = -autoCorr[2]/autoCorr[1];
  coef[1] = refCoef[1];
  Ecopy = (1.0-refCoef[1]*refCoef[1])*Ecopy;
  if ((m-1) > 0){
    for(i = 2; i <= m; i++){
      i1 = i-1;
      sum = 0.0;
      for (j = 1; j <= i1; j++){
	ij1 = i-j+1;
	sum = sum+coef[j]*autoCorr[ij1];
      }	   
      refCoef[i] = -(autoCorr[i+1]+sum)/Ecopy;
      coef[i] = refCoef[i];
      for (j = 1; j <= i; j++)
	b[j] = coef[j];
      for (j = 1; j <= i1; j++) {
	imj = i-j;
	coef[j] = coef[j]+refCoef[i]*b[imj];
      }
      Ecopy = Ecopy*(1.0-refCoef[i]*refCoef[i]);
    }
  }
  Ecopy = Ecopy/autoCorr[1];
  *resEnergy = Ecopy;
}
/****************************************************************************
   RemoveAverage  - removes dc component of a signal             */

 void RemoveAverage(float *derv,int nfft,float *ave)
 {
    int i;
 
    *ave=0.0;
    for(i=1;i<=nfft;i++) *ave+=derv[i];
    *ave = *ave/(float)nfft;
    for(i=1;i<=nfft;i++) derv[i]=derv[i]-*ave;
 }
/****************************************************************************
   ComputeAverage  - Computes the average of a signal             */

 void ComputeAverage(float *signal,int npts,float *ave)
 {
    int i;
    
    *ave=0.0;
    for (i = 1; i <= npts ; i++) {
      *ave = *ave + signal[i];
    }
    *ave = *ave/(float)npts;
 }

/*******************************************************************************
*	LpAnal performs LP analysis on a signal

*	Inputs -
*	signal - array of length npts

*       outputs - 
*	resEnergy - residual signal of length npts
*	coef - array of AR coefficients of length order
*	gain - energy
**************************************************************************** */

  void LpAnal(float *signal,float *resEnergy,int npts, 
	      float *coef,int order,float *gain){

    static float 	*sw, *swOrg, *refCoef, w, ee;
    static double       *autoCorr;
    int 	        i,k,j,lim;
    static int          flag = 0;
    static int          order1;
    if (flag == 0) {
      sw = (float *) AllocFloatArray(sw, npts+order+1);
      swOrg = (float *) AllocFloatArray(swOrg, npts+order+1);
      order1 = order+1;
      autoCorr = (double *) calloc(order1+1, sizeof(double));
      refCoef = (float *) AllocFloatArray(refCoef,order1+1);
      flag = 1;
      for (i = 1; i <= order; i++) {
        sw[i] = 0;
        swOrg[i] = 0;
      }
    }
    for (i = 1; i <= npts; i++) {
      sw[i+order] = signal[i]*HamDw(i,npts);
      swOrg[i+order] = signal[i];
    }

    for (k = 1; k <= order1; k++){
      i = k - 1;
      autoCorr[k] = 0.0;
      lim = npts-i;
      for (j = 1; j <= lim; j++)
	autoCorr[k] = autoCorr[k]+sw[j+order]*sw[j+order+i];
    }
    Durbin(autoCorr,order,coef,refCoef,&ee);
    *gain = sqrt(ee*autoCorr[1]);
    /*        printf("gain = %f\n", *Gain); */
    for (j = 1; j <= npts; j++) {
       resEnergy[j] = swOrg[j+order];
      for (k = 1; k <= order; k++)  
	resEnergy[j] = resEnergy[j] + coef[k]*swOrg[j+order-k];
    }
    for (i = 1; i <= order; i++) {
      sw[i] = sw[npts+i];
      swOrg[i] = swOrg[npts+i];
    }
  }
/*****************************************************************************
 Routine LPSpectrum computes the LP magnitude and phase Spectrum from the
 AR coefficients. input a contains the coefficients, order is the LP order, 
 gain is the residual gain from lp analysis and nfft and mfft are the FFT order
 parameters (nfft = 2**mfft). Outputs are obtained in mag and phase arrays.

 Inputs -
         a[] - ARcoeffs of the ARmodel whose Spectrum is reqd.
         order - LPorder (AR order)
         nfft - DFT order (nfft - LPorder = no. of zeros to be padded) 
         mfft - nfft = 2**mfft
         gain - gain of the ARmodel (constant in the numerator)
 Outputs -
         mag[] - array of DFT magnitude
         phase[] - array of DFT phase

******************************************************************************/
void LPSpectrum(a,order,mag,phase,nfft,mfft,gain)
int              nfft,mfft,order;
float            a[],mag[],phase[],gain;
{
  static float          *ac, *ax, *ay;
  int                   i,nsign;
  static int            flag=0;
  
  if (flag == 0) {
    ac = (float *) AllocFloatArray(ac, nfft+1);
    ax = (float *) AllocFloatArray(ax, nfft+1);
    ay = (float *) AllocFloatArray(ay, nfft+1);
    flag = 1;
  }
  ac[1] = 1.0;
  for (i=1; i <= order; i++) 
    ac[i+1] = a[i];
  
  for (i=order+2;i<=nfft;i++) 
    ac[i] = 0;
  nsign = -1;
  Rfft(ac,ax,ay,mfft,nfft,nsign);
  SpectrumReal(nfft,ax,ay,mag,phase);
  /*  printf("gain = %f\n", gain); */
  fflush(stdout);
  for (i=1;i<=nfft;i++) {
    if (mag[i] != 0)
      mag[i]=(gain/mag[i]);
    else
      mag[i] = 1.0E-10;
     /*printf("magaft %d = %f\n",i,mag[i]); */
  }
}

/*****************************************************************************
 Routine LogSpectrum computes the Log Spectrum 
 of the input spectral coefficients (in place).

 Inputs -
         Spectrum[] - array of spectral coefficients
         npts - number of spectral coefficients) 

Outputs - Spectrum[] - contains the log Spectrum
******************************************************************************/
void LogSpectrum(float Spectrum[],int npts)
{
   int i;
   double temp;

   for(i=1; i<=npts; i++){
      temp = ((double)Spectrum[i]*(double)Spectrum[i]);
      if(temp < 0.00000001){
         temp = 0.00000001;
         puts("\n WARNING from LogSpectrum: Squared magn < 0.00000001");
      }
      Spectrum[i] = (float)(10.0*log10(temp));
   }
}

/************************************************************************
  CepSmooth smoothes magnitude amag to give smthAmag. Order nfft.
              2**mfft = nfft. winlen is Window size.

  Inputs -
          amag[] - array of nfft spectral coefficients
          nfft - FFT order
          mfft - 2**mfft = nfft
	  winlen - cepstral Window length
  Outputs -
          smthAmag[] - smoothed Spectrum of nfft coefficients
	  co - c0 cepstral value
*************************************************************************/
void CepSmooth(float amag[], float smthAmag[], int mfft, int nfft,  
	       int winlen,float *c0, float gamma)
{
  static              int flag = 0;
  static float        *resAx, *resAy, *resAxcop,
                      *resAmag, *resPhase;
  int i;
  if (flag == 0) {
    resAx = (float *) AllocFloatArray(resAx, nfft+1);
    resAy = (float *) AllocFloatArray(resAy, nfft+1);
    resAxcop = (float *) AllocFloatArray(resAxcop, nfft+1);
    resAmag = (float *) AllocFloatArray(resAmag, nfft+1);
    resPhase = (float *) AllocFloatArray(resPhase, nfft+1);
    flag = 1;
  } 
 *c0 = 0;
 for ( i=1; i<=nfft; i++ ) {
   *c0 = *c0 + amag[i]*amag[i];
   if (amag[i] <= 0.000001) 
     resAmag[i] = 0.000001;
   else
     resAmag[i] = exp(gamma*log(amag[i]));
 }
 *c0 = sqrt(*c0/nfft);
 Rfft(resAmag, resAx, resAy, mfft, nfft, 1);

 for ( i=2; i<=winlen; i++ )
    {
     resAxcop[i] = resAx[i]*HanW(i,winlen);
     resAxcop[nfft-i+2] = resAxcop[i];
    }
 resAxcop[1] = resAx[1];
 for(i=winlen+1;i<=nfft-winlen+1;i++) 
   resAxcop[i] = 0.0;
 
 Rfft(resAxcop, resAmag, resAy, mfft, nfft, -1);
 SpectrumReal(nfft, resAmag, resAy, smthAmag, resPhase);
}

/******************************************************************
   ComputeCepstrum computes the cepstrum 
Input :
    amag : magnitude Spectrum
    mfft, nfft : order of FFT
Output :
    cepstrum : an array of length nfft/2
******************************************************************/
void ComputeCepstrum(float amag[],float cepstrum[],int mfft, int nfft)
{
 float resAx[1025],resAy[1025],resAmag[1025];
 int i,nfby2;
 nfby2 = nfft/2;
 for ( i=1; i<=nfft; i++ ) resAmag[i] = log(amag[i]);
 Rfft(resAmag,resAx,resAy,mfft,nfft,1);
 for( i = 1; i < nfby2; i++)
   cepstrum[i] = resAx[i];
}

/*****************************************************************
 *  function GeneratePseudoDct generates the discrete cosine
transform 
Inputs - numRows and numColumns
outputs - VECTOR_OF_F_VECTORS melCosineTransform
*******************************************************************/

VECTOR_OF_F_VECTORS  *GeneratePseudoDct (int offset, int numRows, int numColumns)
{
  int                     i, j;
  VECTOR_OF_F_VECTORS     *vfv;
  float                   *list;

  vfv = (VECTOR_OF_F_VECTORS *) calloc(numRows, sizeof(F_VECTOR));
  if (vfv == NULL) {
    printf("unable to allocate vfv array \n");
    fflush(stdout);
    exit(-1);
  }
  for (i = offset; i <= offset+numRows-1; i++){
    vfv[i-offset] = (F_VECTOR *) malloc(1);
    if (vfv[i-offset] == NULL) {
      printf("unable to allocate vfv FVector\n");
      fflush(stdout);
      exit(-1);
    }
    list = AllocFloatArray(list,numColumns);
    fflush(stdout);
    vfv[i-offset]->numElements = numColumns;
    vfv[i-offset]->array = list;
    for (j = 1; j <= numColumns; j++){
      vfv[i-offset]->array[j-1] = (double) cos (PI * (i - 0.5)*j/numColumns);
    }
  }
  return(vfv);
}
/***********************************************************************
  
  function Warp Warps an input frequency fin

  Inputs -  fin - input frequency
            minFrequency, maxFrequency - range of frequency scale
	    WarpConstant - Warping constant (0 - 1.0, 0 => no Warping
	                                               1 => heavy Warping
   Outputs - the Warped frequency
**********************************************************************/

float Warp (float fin, float minFrequency, float maxFrequency, 
	    float warpConstant)
{ 
  float                warpedFrequency;
  float                range,value;

  range = maxFrequency - minFrequency;
  value = (fin - minFrequency)*PI/range;
  warpedFrequency = fin + range *(2.0/PI)* 
    (atan(-warpConstant*sin(value)
	  /(1 + warpConstant *value)));
  return(warpedFrequency); 
}

/************************************************************************
TrapezoidalFilter outputs the weights for the trapezoidal filter

     Inputs :
     startFreq - startFrequency of filter
     endFreq   - endFrequency of filter
     fin - input frequency
     trapRatio - trapezoidal ratio

***********************************************************************/

 float TrapezoidalFilter(float startFreq, float endFreq, 
			  float fin, float trapRatio) {
   
   float                  limits,slope;
   float                  weight;
   float                  leftLim, rightLim;

   float centreFreq = startFreq + fabs(endFreq - startFreq)/2.0;
   limits = (centreFreq - startFreq)*(1-trapRatio);
   if(limits != 0){
     slope = 1/limits;
     leftLim = startFreq + limits;
     rightLim = endFreq - limits;
     if (fin < leftLim) 
       weight = (fin-startFreq)*slope;
     else if (fin > rightLim)
       weight = (endFreq-fin)*slope;
     else weight = 1;
   }
   else weight = 1;
   return(weight);
 }

/***********************************************************************

GenerateMelFilters generates a set of melFilters

      Inputs : asdf - the structure containing all the relevant
               input information
      Outputs : filterbank is stored in asdf itself

**********************************************************************/

void GenerateMelFilters (ASDF *asdf)
{
 VECTOR_OF_F_VECTORS           *filterbankWeights;
 int                           *dftIndices;
 float                         frequencyInc;
 float                         dftLineFrequency;
 int                           minDftFrequency,maxDftFrequency;
 float                         *startFrequency,*endFrequency,
                               *centreFrequency;
 int                           *startIndex,*endIndex,*centreIndex;
 int                           i,j, numElem;

 dftIndices = (int *) calloc(asdf->numFilters, sizeof(int));
 if (dftIndices == NULL) {
   /*     printf("problems\n"); */
   fflush(stdout);
   exit(-1);
 }
 frequencyInc = (asdf->maxFrequency - asdf->minFrequency)/(asdf->numFilters+1);
 dftLineFrequency = (int) asdf->samplingRate/asdf->fftSize;
 minDftFrequency = (int) asdf->minFrequency/dftLineFrequency;
 maxDftFrequency = (int) asdf->maxFrequency/dftLineFrequency;
 startFrequency = (float *) calloc(asdf->numFilters, sizeof(float));
 endFrequency = (float *) calloc(asdf->numFilters, sizeof(float));
 centreFrequency = (float *) calloc(asdf->numFilters, sizeof(float));
 
 startIndex = (int *) calloc(asdf->numFilters, sizeof(int));
 endIndex = (int *) calloc(asdf->numFilters, sizeof(int));
 centreIndex = (int *) calloc(asdf->numFilters, sizeof(int));
 for (i = 0; i < asdf->numFilters; i++) 
   centreFrequency[i] = asdf->minFrequency + (i+1)*frequencyInc;
 
 if (asdf->bandwidthScale == 0.0)
   for (i = 0; i < asdf->numFilters; i++) {
     startFrequency[i] =centreFrequency[i] - frequencyInc;
     endFrequency[i] = centreFrequency[i] + frequencyInc;
     /*      printf("sf = %f cf = %f ef = %f\n", startFrequency[i], centreFrequency[i], endFrequency[i]); */
   } else
     for (i = 0; i < asdf->numFilters; i++) {
       if (centreFrequency[i] < 1000) {
	 startFrequency[i] = centreFrequency[i]- asdf->bandwidthScale*137.50;
	 endFrequency[i] = centreFrequency[i] + asdf->bandwidthScale*137.50;
       } else if ((centreFrequency[i] >= 1000) && 
		  (centreFrequency[i] <= 0.4*asdf->samplingRate)){
	 startFrequency[i] = centreFrequency[i] - asdf->bandwidthScale*1.11*
	   (exp(log(centreFrequency[i])*0.7));
	 endFrequency[i] = centreFrequency[i] + asdf->bandwidthScale*1.11*
	   (exp(log(centreFrequency[i])*0.7));
       } else {
	 startFrequency[i] = centreFrequency[i] - asdf->bandwidthScale*10.84*
	   (exp(log(centreFrequency[i])*0.4));
	 endFrequency[i] = centreFrequency[i] + asdf->bandwidthScale*10.84*
	   (exp(log(centreFrequency[i])*0.4));
       }
       /*      printf("sf = %f cf = %f ef = %f\n", startFrequency[i], centreFrequency[i], endFrequency[i]); */
     }
 for (i = 0; i < asdf->numFilters; i++) {
   startIndex[i] = (int) (Warp(startFrequency[i],asdf->minFrequency,
				asdf->maxFrequency,
			       asdf->filterWarp)/dftLineFrequency);
   endIndex[i] = (int ) (Warp(endFrequency[i],asdf->minFrequency,
			       asdf->maxFrequency,
			      asdf->filterWarp)/dftLineFrequency);
   centreIndex[i] = (int)(Warp(centreFrequency[i], asdf->minFrequency,
				asdf->maxFrequency,
			       asdf->filterWarp)/dftLineFrequency);
 }
 filterbankWeights = (VECTOR_OF_F_VECTORS *) 
   malloc(asdf->numFilters*sizeof(F_VECTOR *));
 if (filterbankWeights == NULL) {
   printf("problems allocating filterbankWeights");
   exit(-1);
 }
 for (i = 0; i < asdf->numFilters; i++) {
   numElem = endIndex[i] - startIndex[i]+1;
   filterbankWeights[i] = (F_VECTOR *) AllocFVector(numElem);
   dftIndices[i] = startIndex[i];
   for (j = startIndex[i]; j<= endIndex[i]; j++){
     filterbankWeights[i]->array[j-startIndex[i]] = 
       TrapezoidalFilter(startIndex[i], endIndex[i], j, 
			 asdf->trapezoidalRatio);
     /*     printf("fbWts %d %d = %f\n",i, 
	    j-startIndex[i],filterbankWeights[i]->array[j-startIndex[i]]); */
   }
   
 }
 asdf->filterbankWeights = filterbankWeights;
 asdf->dftIndices = dftIndices;
}


/*-------------------------------------------------------------------------
 *  LinearTransformationOfFVector performs  a DCT on a F_VECTOR and
 *  outputs an F_VECTOR
 
    Inputs :
         inVect : a vector of type F_VECTOR
         melCepstrumCosineTransform : DCT
    Outputs :
         outVect : a vector of type F_VECTOR
---------------------------------------------------------------------------*/
    
void LinearTransformationOfFVector (F_VECTOR *inVect, VECTOR_OF_F_VECTORS 
				       *melCepstrumCosineTransform, 
				       int numRows, int numCols,
				       F_VECTOR *outVect) {

int                                    i, j;

for (i = 0; i < numRows; i++){
  outVect->array[i] = 0;
  for (j = 0; j < numCols; j++)
    outVect->array[i] = outVect->array[i] + 
      inVect->array[j]*melCepstrumCosineTransform[i]->array[j];
}
}

/*----------------------------------------------------------------------------

  FilterbankEnergyIntegration : integrates the given Spectrum with a comb
filter.

  Inputs :
           asdf - structure which contains the front-end-parameters
           Spectrum - Spectrum of the given frame of speech

  Outputs :
           fvect - contains the output filterbank energies.

-----------------------------------------------------------------------------*/
F_VECTOR *FilterbankEnergyIntegration(ASDF *asdf, float *Spectrum,
				      F_VECTOR *fvect) {

  static int                 mfft,nfft,WindowSize,numFilters,nfby2;
  static ASDF                *prevAsdf=NULL;
  static VECTOR_OF_F_VECTORS *filterbankWeights;
  static int                 *dftIndices;
  float                      sum;
  int                        i, j;

  if (prevAsdf != asdf) {
    mfft = GetIAttribute(asdf, "fftOrder");
    nfft = GetIAttribute(asdf,"fftSize");
    nfby2 = nfft/2;
    numFilters = GetIAttribute(asdf,"numFilters");
    filterbankWeights = asdf->filterbankWeights;
    dftIndices = asdf->dftIndices;
    prevAsdf = asdf;
  }
  
  fvect->numElements = numFilters;
  for (i = 0; i < numFilters; i++) {
    fvect->array[i] = 0;
    for (j = dftIndices[i]; 
	 j < dftIndices[i]+filterbankWeights[i]->numElements; j++) {
      if (j <= nfby2){
        fvect->array[i] = fvect->array[i] + Spectrum[j]*Spectrum[j]* 
        filterbankWeights[i]->array[j-dftIndices[i]];
	/*          printf("fbWts %d %d  = %f spec %d = %f dftInd %d = %d\n",
     i, j-dftIndices[i], filterbankWeights[i]->array[j-dftIndices[i]], 
      j,Spectrum[j], i, dftIndices[i]);*/ 
      }
    }
  }
  return(fvect);
}
/*--------------------------------------------------------------------------

 gets a frame frameIndex of speech from the waveform waveform 

---------------------------------------------------------------------------*/


float *FrameCompWaveform(short *waveform, float *array, int frameIndex, 
			 int frameShift, int frameSize, long samples) {

int           i;

  for (i = 0; i < frameSize; i++)
    if (frameShift*frameIndex+i < samples) 
      array[i+1] = waveform[frameShift*frameIndex+i];
    else
      array[i+1] = 0.0;
  return(array);
}   



/*---------------------------------------------------------------------------

Computes the energy in an array of size frameSize

Inputs : a list array of frameSize

Outputs : energy

-----------------------------------------------------------------------------*/


float ComputeEnergy (float *array, int frameSize) {
int i;
float sum = 0;
for (i = 1; i <= frameSize; i++){
  sum = sum + array[i]*array[i];
}
sum = sum/frameSize;
return(sum);
}

/*-----------------------------------------------------------------------------

Computes the zero-crossing rate in array of size of frameSize

Inputs : a list array of frameSize

Outputs : zero-crossing rate
------------------------------------------------------------------------------*/

float ComputeZeroCrossing(float *array, int frameSize) {
  float              *signum;
  int                i;
  float              sum;

  signum = (float*) calloc(frameSize+1, sizeof(float));
  for (i = 1; i <= frameSize; i++)
    if (array[i] > 0)
      signum[i] = 1;
    else
      signum[i] = -1;
  for(i = 1; i < frameSize; i++)
    sum = sum + fabs(signum[i] - signum[i-1]);
  sum = sum/(2*frameSize);
return(sum);
}

/*-----------------------------------------------------------------------------

Estimates the spectral flatness in an array of size frameSize 

Inputs : a list array of frameSize, another list residual of frameSize

Outputs : spectral flatness

------------------------------------------------------------------------------*/
float ComputeSpectralFlatness(float *array, float *residual, int frameSize){ 
  float                  sigEnergy, resEnergy;
  float                  coef[11];
  float                  gain;
  float                  flatness;

  sigEnergy = ComputeEnergy(array, frameSize);
  LpAnal(array,residual,frameSize, coef, 10,&gain);
  resEnergy = ComputeEnergy(residual, frameSize);
  flatness = resEnergy/sigEnergy;
  return(flatness);
}

/*------------------------------------------------------------------------------

 Classifies a given speech signal into voiced and unvoiced frames 

Input :

waveform      - given waveform
samples       - number of samples
frameShift   - shift in location of sliding Window
frameSize    - size of a frame of speech

Output : 
vU           - a binary array 
              - a 1 in a particular location indicates that the given frame of
                speech is voiced
              - a 0 in a particular location indicates that the given frame of 
                speech is unvoiced

-----------------------------------------------------------------------------*/


void VoicedUnvoiced(short *waveform, long samples, short *vU, 
		    int frameShift, int frameSize) {
  float             *zeroCrossings, *energy, *specFlatness, *residual;
  int               numFrames,i;
  float             *array;
  float             medZero, medEnergy, medSpecFlatness;
  float             ThZero = 0.5, ThEnergy = 0.5, ThSpec = 0.5;

  numFrames = (int) samples/frameShift;
  zeroCrossings = (float *) calloc(numFrames+1,sizeof(float));
  energy = (float *) calloc(numFrames+1,sizeof(float));
  specFlatness = (float *) calloc(numFrames+1, sizeof(float)); 
  array = (float *) calloc(frameSize+1,sizeof(float));
  residual = (float *) calloc(frameSize+1,sizeof(float));
  for (i = 1; i <= numFrames; i++){
    array = FrameCompWaveform(waveform, array, i-1, 
			      frameShift, frameSize, samples);
    
    zeroCrossings[i] = ComputeZeroCrossing(array, frameSize); 
    energy[i] = ComputeEnergy(array,frameSize);
    specFlatness[i] = ComputeSpectralFlatness(array,residual, frameSize);
    /*  printf("zero %f energy %f flat %f \n", zeroCrossings[i], energy[i], specFlatness[i]);
	scanf("*c"); */
  }
  medZero = Median(zeroCrossings,numFrames);
  medEnergy = Median(energy,numFrames);
  medSpecFlatness = Median(specFlatness, numFrames);
  /*  printf("numFrames %d maxzero %f maxenergy %f maxflat %f \n",numFrames, medZero, medEnergy,medSpecFlatness);
      printf("numFrames %d thzero %f thenergy %f thflat %f \n",numFrames, medZero, medEnergy, medSpecFlatness); */
  for (i = 1; i <= numFrames; i++)
    if ((zeroCrossings[i] > medZero) &&
	(energy[i] < medEnergy) &&
	(specFlatness[i] > medSpecFlatness))
      vU[i] = 0;
    else
      vU[i] = 1;
}


/****************************************************************************
*	The following subroutine computes the Group delay Cepstrum
*       of a given signal
*
*	inputs : signal npts 
*	mfft : fft order nfft = 2**mfft
*	winlen : number of cepstral coeffts
*
*	Outputs :
*	groupDelay(nfft) - Group delay function
*****************************************************************************
*
*/
float *GroupDelayCepstrum(float *signal,int npts, int nfft, int mfft, 
			  int winlen, float *groupDelayCepstrum){
  static float		*sigTemp;
  static float		*ax, *ay,*amag,*phase, *groupDelay;
  static int            nfBy2,flag = 0;
  static complex 	*cSig,*cfSig;
  static complex	*cfTemp1,*cfTemp2,u;
  static float          epsilon = 0.0001;
  int		        i;
  float                 Ave, c0;

  if (flag == 0) {
    nfBy2 = nfft/2;
    ax = (float *) AllocFloatArray(ax,nfft+1);
    ay = (float *) AllocFloatArray(ay,nfft+1);
    amag = (float *) AllocFloatArray(amag,nfft+1);
    phase = (float *) AllocFloatArray(phase,nfft+1);
    groupDelay = (float *) AllocFloatArray(groupDelay,nfft+1);
    sigTemp = (float *) AllocFloatArray(sigTemp,nfft+1);
    cSig = (complex *) calloc (nfft+1, sizeof(complex));
    cfSig = (complex *) calloc (nfft+1, sizeof(complex));
    cfTemp1 = (complex *) calloc (nfft+1, sizeof(complex));
    cfTemp2 = (complex *) calloc (nfft+1, sizeof(complex));
    if ((cSig == NULL) || (cfSig == NULL) || (cfTemp1 == NULL) ||
	(cfTemp2 == NULL)) {
      printf("unable to allocate complex array in Group Delay\n");
      fflush(stdout);
      exit(-1);
    } 
    flag = 1;
  }
  for (i = 1; i <= npts; i++)
    sigTemp[i] = signal[i];
  for(i = npts+1; i<= nfft; i++)
    sigTemp[i] = 0.0;
  for (i = 1; i <= npts; i++){
    u.re = sigTemp[i];
    u.im = 0.0;
    cSig[i] = u;
    u.re = (float)(i-1);
    u.im = 0.0;
    cmul(cfSig[i],u,cSig[i]);
  }
  for (i = npts+1; i <= nfft; i++) {
    cSig[i].re = 0.0;
    cSig[i].im = 0.0;
    cfSig[i].re = 0.0;
    cfSig[i].im = 0.0;
  }
  Rfft(sigTemp,ax,ay,mfft,nfft,-1);
  SpectrumReal(nfft,ax,ay,amag,phase);
  c0 = 0;
  for (i = 1; i <= nfft; i++)
    c0 = c0 + amag[i]*amag[i];
  c0 = sqrt(c0/nfft);
  Cfft(cSig,cfTemp1,mfft,nfft,-1);
  Cfft(cfSig,cfTemp2,mfft,nfft,-1);
  for (i = 1; i <= nfft; i++){
    conjg(u,cfTemp1[i]);
    cmul(cfTemp2[i],cfTemp2[i],u);
    u.re = cfTemp2[i].re;
    u.im = cfTemp2[i].im;
    if (amag[i] > epsilon) {
      cfTemp2[i].re = u.re/(amag[i]*amag[i]);
      cfTemp2[i].im = u.im/(amag[i]*amag[i]);
    } else 
      cfTemp2[i] = cfTemp2[i-1];
    groupDelay[i] = cfTemp2[i].re;
  }
  RemoveAverage(groupDelay,nfft,&Ave); 
  Rfft(groupDelay,ax,ay,mfft,nfft,1);
  for (i = 1; i <= winlen; i++) 
    groupDelayCepstrum[i-1] = ax[i+1]*HanW(i,winlen+1);
  return(groupDelayCepstrum);
}
/****************************************************************************
*	The following subroutine computes the smoothed group delay function of signal
*
*	inputs : signal(npts) 
*	mfft : fft order  = 2**mfft
*	winlen : window length for smoothing
*
*	Outputs :
*	groupDelay(nfft) - Group delay function
*****************************************************************************
*
*/
float *GroupDelay(float *signal,int npts, int nfft, int mfft, int winlen,
		float *groupDelay){
  static float		*ax, *ay,*amag,*phase, *gdCepstrum;
  static int            nfBy2,flag = 0;
  int		        i;
  float                 Ave;

  if (flag == 0) {
    nfBy2 = nfft/2;
    ax = (float *) AllocFloatArray(ax,nfft+1);
    ay = (float *) AllocFloatArray(ay,nfft+1);
    amag = (float *) AllocFloatArray(amag,nfft+1);
    phase = (float *) AllocFloatArray(phase,nfft+1);
    gdCepstrum = (float *) AllocFloatArray(gdCepstrum,winlen);
    flag = 1;
  }
  gdCepstrum = (float *) GroupDelayCepstrum(signal, npts, nfft, mfft,
					    winlen, gdCepstrum);
  ax[1] = 0.0;
  for (i = 2; i <= winlen; i++) {
    ax[i] = gdCepstrum[i-2];
    ax[nfft-i+2] = ax[i];
  }
  for (i = winlen+1; i <= nfft-(winlen-1); i++)
    ax[i] = 0.0;
  Rfft(ax,groupDelay,ay,mfft,nfft,-1);
  return(groupDelay);
}

/****************************************************************************
*	The following subroutine computes the standard group delay function of signal
*
*	inputs : DATA npts points long
*	mfft : fft stages nfft = 2**mfft
*	winlen : window length for zero spectrum
*
*	Outputs :
*	groupDelay(nfft) - Group delay function
*****************************************************************************
*
*/
	void StdGroupDelay(float *signal,int npts, int nfft, int mfft, 
		    float *groupDelay){
	static float		*sigTemp;
	static float		*ax, *ay,*amag,*phase;
        static int              nfBy2,flag = 0;
	static complex 	        *cSig,*cfSig;
	static complex		*cfTemp1,*cfTemp2,u;
	int		        i;

        if (flag == 0) {
	  nfBy2 = nfft/2;
          ax = (float *) calloc(nfft+1, sizeof(float));
          ay = (float *) calloc(nfft+1, sizeof(float));
          amag = (float *) calloc(nfft+1, sizeof(float));
          phase = (float *) calloc(nfft+1, sizeof(float));
          cSig = (complex *) calloc (nfft+1, sizeof(complex));
          cfSig = (complex *) calloc (nfft+1, sizeof(complex));
          cfTemp1 = (complex *) calloc (nfft+1, sizeof(complex));
          cfTemp2 = (complex *) calloc (nfft+1, sizeof(complex));
          sigTemp = (float *) calloc(nfft+1, sizeof(float));
          flag = 1;
	}
	for (i = 1; i <= npts; i++)
	  sigTemp[i] = signal[i];
	for(i = npts+1; i<= nfft; i++)
	  sigTemp[i] = 0.0;
	for (i = 1; i < npts; i++){
	  u.re = sigTemp[i];
          u.im = 0.0;
	  cSig[i] = u;
	  u.re = (float)(i-1);
          u.im = 0.0;
	  cmul(cfSig[i],u,cSig[i]);
	}
	for (i = npts+1; i <= nfft; i++) {
	  cSig[i].re = 0.0;
          cSig[i].im = 0.0;
	  cfSig[i].re = 0.0;
          cfSig[i].im = 0.0;
	}
	Rfft(sigTemp,ax,ay,mfft,nfft,-1);
	SpectrumReal(nfft,ax,ay,amag,phase);
	Cfft(cSig,cfTemp1,mfft,nfft,-1);
	Cfft(cfSig,cfTemp2,mfft,nfft,-1);
	for (i = 2; i <= nfBy2; i++){
          conjg(u,cfTemp1[i]);
	  cmul(cfTemp2[i],cfTemp2[i],u);
	  u.re = cfTemp2[i].re;
          u.im = cfTemp2[i].im;
          cfTemp2[i].re = (u.re/(amag[i]*amag[i]));
          cfTemp2[i].im = (u.im/(amag[i]*amag[i]));
	  groupDelay[i] = cfTemp2[i].re;
	}
        conjg(u,cfTemp1[1]);
	cmul(cfTemp2[1],cfTemp2[1],u);
        if (cabs2(cfTemp1[1]) == 0.0){ 
	  cfTemp2[1].re = cfTemp2[1].re/0.0001;
          cfTemp2[1].im = cfTemp2[1].im/0.0001;
	}
	else{
          cfTemp2[1].re = cfTemp2[1].re/(amag[1]*amag[1]);
          cfTemp2[1].im = cfTemp2[1].im/(amag[1]*amag[1]);
	}
	groupDelay[1] = cfTemp2[1].re;
	}

/****************************************************************************
*	The following subroutine computes the group delay function of 
*       LP Residual of a given signal.
*	inputs : DATA npts points long
*	mfft : fft stages nfft = 2**mfft
*	winlen : window length for zero spectrum
*
*	Outputs :
*	groupDelay(nfft) - Group delay function
*****************************************************************************
*
*/
void LPResGroupDelay(float *signal,int npts, int nfft, int mfft,  
		     int order, float *groupDelay){
  static float		*resEnergy, *coef;
  static float            gain;
  static int              flag = 0;
	
  if (flag == 0) {
    resEnergy = (float *) AllocFloatArray(resEnergy, npts+1);
    coef = (float *) AllocFloatArray(coef, order+2);
    flag = 1;
  }
  LpAnal(signal,resEnergy,npts, coef,order, &gain);
  StdGroupDelay(resEnergy,npts, nfft, mfft, groupDelay);
}


/****************************************************************************
*	The following function computes the minimum phase 
*	group delay cepstrum for the given signal 
*
*	inputs : signal(npts)
*	mfft : fft order  nfft = 2**mfft
*	winlen : window length for minimum phase signal
*	gamma : time compression of minimum phase signal
*
*       Outputs:
*       MinGdCepstrum(winlen) - float *
*       Comment - mean removal not done as MinGd already is zero mean
*****************************************************************************
*
*/
float *MinGdCepstrum(float *signal,int npts, int nfft, int mfft, 
		int winlen, float *minGdCepstrum, float gamma){
  static float	  *sigTemp;
  static int      flag = 0;
  static float	  *ax, *ay,*amag,*phase, *minGd;
  static int      nfBy2;
  int		  i,npeaks;
  
  if (flag == 0) {
    nfBy2 = nfft/2;
    sigTemp = (float *) AllocFloatArray(sigTemp,nfft+1);
    ax =  (float *) AllocFloatArray(ax,nfft+1);
    ay = (float *) AllocFloatArray(ay,nfft+1);
    amag = (float *) AllocFloatArray(amag,nfft+1);
    phase = (float *) AllocFloatArray(phase,nfft+1);
    minGd = (float *) AllocFloatArray(minGd,nfft+1);
    flag = 1;
  }
  for (i = 1; i <= npts; i++)
    sigTemp[i] = signal[i];
  for(i = npts+1; i<= nfft; i++)
    sigTemp[i] = 0.0;
  Rfft(sigTemp,ax,ay,mfft,nfft,-1);
  SpectrumReal(nfft,ax,ay,amag,phase);
  for (i = 1; i <= nfft; i++){
    amag[i] = exp(log(amag[i])*gamma);
  }
  Rfft(amag,sigTemp,ay,mfft,nfft,1);
  for (i = 1; i <= winlen; i++)
    minGdCepstrum[i-1] = sigTemp[i]*HanW(i,winlen+1);
  return(minGdCepstrum);
}

/****************************************************************************
*	The following function computes the Minimum phase 
*	group delay function from the given signal 
*
*	inputs : signal(npts)
*	mfft : fft order  nfft = 2**mfft
*	winlen : window length for minimum phase signal
*	alfa : compression required for magnitude spectrum
*
*       Outputs:
*       MinGd - float *
*****************************************************************************
*
*/
float *MinGd(float *signal,int npts, int nfft, int mfft, 
		int winlen, float *minGd, float gamma){

  static float	  *ax, *ay,*amag,*phase; 
  static float    *sigTemp, *minGdCepstrum;
  static int      flag = 0;
  static int      nfBy2;
  int		  i;
  float           Ave;
  
  if (flag == 0) {
    nfBy2 = nfft/2;
    sigTemp = (float *) AllocFloatArray(sigTemp,nfft+1);
    ax =  (float *) AllocFloatArray(ax,nfft+1);
    ay = (float *) AllocFloatArray(ay,nfft+1);
    amag = (float *) AllocFloatArray(amag,nfft+1);
    phase = (float *) AllocFloatArray(phase,nfft+1);
    minGdCepstrum = (float *) AllocFloatArray(minGdCepstrum, winlen);
    flag = 1;
  }
  minGdCepstrum = MinGdCepstrum(signal, npts, nfft, mfft,
				winlen, minGdCepstrum, gamma);
  for (i = 1; i <= winlen; i++) 
    sigTemp[i] = minGdCepstrum[i-1];
  for (i = winlen+1; i <= nfft; i++)
    sigTemp[i] = 0.0;
  Rfft(sigTemp,ax,ay,mfft,nfft,-1);
  SpectrumReal(nfft,ax,ay,amag,phase);
  for(i = 1; i <=  nfBy2-1; i++) {
    minGd[i] = phase[i] - phase[i+1];
    minGd[nfft-i+1] = minGd[i];
  }
  minGd[nfBy2] = minGd[nfBy2-1];
  minGd[nfBy2+1] = minGd[nfBy2];
  return(minGd);
}

/****************************************************************************
*	The following function computes  the
*	smoothed magnitude spectrum of the given signal
*       Cepstral smoothing is used.
*
*	inputs : signal - npts
*	mfft : fft order nfft = 2**mfft
*	winlen : window length for zero spectrum
*
*	Outputs :
*	modGdCepstrum(winlen) - modified group delay cepstrum
*****************************************************************************
*/
float *SmoothMagSpectrum(float *signal, int npts, int nfft, int mfft, int smthWinSize,
                         float *smthMag) {
  static float		*ax, *ay, *amag, *phase, *modGd, *sigTemp;
  static float		*cepAmag, c0;
  static int            nfBy2, flag = 0, i;
  if (flag == 0) {
    nfBy2 = nfft/2;
    ax = (float *) calloc(nfft+1, sizeof(float));
    ay = (float *) calloc(nfft+1, sizeof(float));
    amag = (float *) calloc(nfft+1, sizeof(float));
    cepAmag = (float *) calloc(nfft+1, sizeof(float));
    phase = (float *) calloc(nfft+1, sizeof(float));
    sigTemp = (float *) calloc(nfft+1, sizeof(float));
    flag = 1;
  }

 
  for (i = 1; i <= npts; i++)
    sigTemp[i] = signal[i];
  for(i = npts+1; i<= nfft; i++)
    sigTemp[i] = 0.0;
 Rfft(sigTemp,ax,ay,mfft,nfft,-1);
  SpectrumReal(nfft,ax,ay,amag,phase);
  for( i = 1; i <= nfft; i++)
    ax[i] = amag[i];
  CepSmooth(ax,cepAmag,mfft,nfft, smthWinSize, &c0, 1.0);
  for (i = 0; i <= nfBy2; i++)
    smthMag[i] = cepAmag[i+1];
  return(smthMag);
  }
/****************************************************************************
*	The following function computes  the
*	modified group delay cepstrum (ncn) for the given signal
*
*	inputs : signal - npts
*	mfft : fft order nfft = 2**mfft
*	winlen : window length for zero spectrum
*
*	Outputs :
*	modGdCepstrum(winlen) - modified group delay cepstrum
*****************************************************************************
*/
float *ModGdCepstrumNcN(float *signal,int npts,int nfft,int mfft, 
		     int winlen, int smthWinSize, float *modGdCepstrum, 
		     float alfaP, float alfaN,  
		     float gamma, int gdsign, 
		     int removeLPhase, int removeMin,  
                     int startIndex, int endIndex) { 
  static float		*sigTemp;
  static float		*ax, *ay, *amag, *phase, *modGd;
  static float		*cepAmag;
  static int            nfBy2, flag = 0;
  static complex        *cSig, *cfSig;
  static complex	*cfTemp1, *cfTemp2, u;
  int		        i, sign;
  float		        c0, Ave, rmax, abscfTemp, logcfTemp, alfaLogTemp;
  float                 minValue;


  if (flag == 0) {
    nfBy2 = nfft/2;
    ax = (float *) calloc(nfft+1, sizeof(float));
    ay = (float *) calloc(nfft+1, sizeof(float));
    amag = (float *) calloc(nfft+1, sizeof(float));
    cepAmag = (float *) calloc(nfft+1, sizeof(float));
    phase = (float *) calloc(nfft+1, sizeof(float));
    modGd = (float *) calloc(nfft+1, sizeof(float));
    cSig = (complex *) calloc (nfft+1, sizeof(complex));
    cfSig = (complex *) calloc (nfft+1, sizeof(complex));
    cfTemp1 = (complex *) calloc (nfft+1, sizeof(complex));
    cfTemp2 = (complex *) calloc (nfft+1, sizeof(complex));
    sigTemp = (float *) calloc(nfft+1, sizeof(float));
    flag = 1;
  }
  for (i = 1; i <= npts; i++)
    sigTemp[i] = signal[i];
  for(i = npts+1; i<= nfft; i++)
    sigTemp[i] = 0.0;
  for (i = 1; i < npts; i++){
    u.re = sigTemp[i];
    u.im = 0.0;
    cSig[i] = u;
    u.re = (float)(i-1);
    u.im = 0.0;
    cmul(cfSig[i],u,cSig[i]);
  }
  for (i = npts+1; i <= nfft; i++) {
    cSig[i].re = 0.0;
    cSig[i].im = 0.0;
    cfSig[i].re = 0.0;
    cfSig[i].im = 0.0;
  }
  Cfft(cSig,cfTemp1,mfft,nfft,-1);
  Cfft(cfSig,cfTemp2,mfft,nfft,-1);
  Rfft(sigTemp,ax,ay,mfft,nfft,-1);
  SpectrumReal(nfft,ax,ay,amag,phase);
  for( i = 1; i <= nfft; i++)
    ax[i] = amag[i];
  CepSmooth(ax,cepAmag,mfft,nfft, smthWinSize, &c0, gamma);
  for (i = 1; i <= nfBy2+1; i++){
    conjg(u,cfTemp1[i]);
    cmul(cfTemp2[i],cfTemp2[i],u);
    u.re = cfTemp2[i].re;
    u.im = cfTemp2[i].im;
    cfTemp2[i].re = (u.re/(cepAmag[i]*cepAmag[i]));
    cfTemp2[i].im = (u.im/(cepAmag[i]*cepAmag[i]));
    modGd[i] = cfTemp2[i].re;
    if (i > 1)
      modGd[nfft-i+2] = modGd[i];
  }
  /* Remove DC component of modGd */

  /* RemoveAverage(modGd,nfft,&Ave); */

  /* Reduce dynamic range of modGd if necessary */
    for (i = 1; i <= nfft; i++) {
      if (modGd[i] < 0) 
	sign = -1;
      else
	sign = 1;
      abscfTemp = fabs(modGd[i]);
      if (abscfTemp == 0)
	modGd[i] = 0;
      else {
	logcfTemp = log(abscfTemp);
	if ((alfaN == 0) && (alfaP == 0))
	  if (gdsign == 1)
	    modGd[i] = sign*logcfTemp;
	  else
	    modGd[i] = logcfTemp;
	else {
	  if (sign > 0.0)
	    modGd[i] = exp(alfaP*logcfTemp);
	  else 
	    modGd[i] = exp(alfaN*logcfTemp);
	  if (gdsign == 1)
	    modGd[i] = sign*modGd[i];
	}
      }
    }  
    /* Bandlimit modGd */
    if (removeLPhase) {
      for (i = 2; i < nfBy2; i++) {
	modGd[i] = modGd[i] - modGd[1];
	modGd[nfft - i + 2] = modGd[i];
      }
      modGd[1] = 0;
    }
    minValue = modGd[IminActual(modGd, nfft)];
    if (removeMin) 
      for (i = 1; i <nfft; i++) 
        modGd[i] = modGd[i] - minValue;
      
    if ((startIndex > 1) || (endIndex < nfBy2)) {
      for (i = 1; i < startIndex; i++) { 
	modGd[i] = 0.0;
	if (i > 1)
	  modGd[nfft - i+2] = modGd[i];
      }
      for (i = endIndex+1; i <= nfBy2; i++) {
	modGd[i] = 0;
	if (i > 1)
	  modGd[nfft - i+ 2] = 0;
      }
    }

  Rfft(modGd,ax,ay,mfft,nfft,1);
  for (i = 1; i <= winlen; i++) 
    modGdCepstrum[i-1] = ax[i];
  return(modGdCepstrum);  
}

/****************************************************************************
*	The following function computes  the
*	modified group delay cepstrum (ncn) for the given signal
*
*	inputs : signal - npts
*	mfft : fft order nfft = 2**mfft
*	winlen : window length for zero spectrum
*
*	Outputs :
*	modGdCepstrum(winlen) - modified group delay cepstrum
*****************************************************************************
*
*/
float *ModGdCepstrum(float *signal,int npts,int nfft,int mfft, 
		     int winlen, int smthWinSize, float *modGdCepstrum, float alfaP, float alfaN,  
		     float gamma, int gdsign, int removeLPhase, int removeMin, int startIndex, int endIndex) { 

static  float                     *modGdCepTemp;
static  int                       flag = 0;
int                               i;
if (flag == 0) {
  modGdCepTemp = (float *) calloc(winlen+1, sizeof(float));
  flag = 1;
}

modGdCepTemp = ModGdCepstrumNcN (signal, npts, nfft, mfft, winlen+1, smthWinSize,
				 modGdCepTemp, alfaP, alfaN, gamma, gdsign, 
                                 removeLPhase, removeMin,
				 startIndex, endIndex);
for (i = 0; i < winlen; i++)
  modGdCepstrum[i] = modGdCepTemp[i+1]/(float)(i+1);
return(modGdCepstrum);
}
/****************************************************************************
*	The following function computes  the
*	modified group delay cepstrum  for the given signal using DCT
*
*	inputs : signal - npts
*	mfft : fft order nfft = 2**mfft
*	winlen : window length for zero spectrum
*
*	Outputs :
*	modGdCepstrum(winlen) - modified group delay cepstrum
*****************************************************************************
*
*/
float *ModGdCepstrum_DCT(float *signal,int npts,int nfft,int mfft, 
		     int winlen, int smthWinSize, float *modGdCepstrum, float alfaP, float alfaN,  
		     float gamma, int gdsign, int removeLPhase, int startIndex, int endIndex, float scaleDCT) { 
  static float		     *sigTemp;
  static float		     *ax, *ay, *amag, *phase, *modGd;
  static float		     *cepAmag;
  static int                 nfBy2, flag = 0;
  static complex             *cSig, *cfSig;
  static complex	     *cfTemp1, *cfTemp2, u;
  static VECTOR_OF_F_VECTORS *discreteCosineTransform;
  static F_VECTOR            *fvect, *modGdFvect;
  int		              i, sign;
  float		              c0, Ave, rmax, abscfTemp, 
                              logcfTemp, alfaLogTemp;
  
  if (flag == 0) {
    nfBy2 = nfft/2;
    ax = (float *) calloc(nfft+1, sizeof(float));
    ay = (float *) calloc(nfft+1, sizeof(float));
    amag = (float *) calloc(nfft+1, sizeof(float));
    cepAmag = (float *) calloc(nfft+1, sizeof(float));
    phase = (float *) calloc(nfft+1, sizeof(float));
    modGd = (float *) calloc(nfft+1, sizeof(float));
    cSig = (complex *) calloc (nfft+1, sizeof(complex));
    cfSig = (complex *) calloc (nfft+1, sizeof(complex));
    cfTemp1 = (complex *) calloc (nfft+1, sizeof(complex));
    cfTemp2 = (complex *) calloc (nfft+1, sizeof(complex));
    sigTemp = (float *) calloc(nfft+1, sizeof(float));
    discreteCosineTransform = (VECTOR_OF_F_VECTORS *) GeneratePseudoDct(1, winlen, nfBy2);
    fvect = (F_VECTOR *) AllocFVector(nfBy2);
    modGdFvect = AllocFVector(winlen);
    flag = 1;
  }
  for (i = 1; i <= npts; i++)
    sigTemp[i] = signal[i];
  for(i = npts+1; i<= nfft; i++)
    sigTemp[i] = 0.0;
  for (i = 1; i < npts; i++){
    u.re = sigTemp[i];
    u.im = 0.0;
    cSig[i] = u;
    u.re = (float)(i-1);
    u.im = 0.0;
    cmul(cfSig[i],u,cSig[i]);
  }
  for (i = npts+1; i <= nfft; i++) {
    cSig[i].re = 0.0;
    cSig[i].im = 0.0;
    cfSig[i].re = 0.0;
    cfSig[i].im = 0.0;
  }
  Cfft(cSig,cfTemp1,mfft,nfft,-1);
  Cfft(cfSig,cfTemp2,mfft,nfft,-1);
  Rfft(sigTemp,ax,ay,mfft,nfft,-1);
  SpectrumReal(nfft,ax,ay,amag,phase);
  for( i = 1; i <= nfft; i++)
    ax[i] = amag[i];
  CepSmooth(ax,cepAmag,mfft,nfft, smthWinSize, &c0, gamma);
  for (i = 1; i <= nfBy2+1; i++){
    conjg(u,cfTemp1[i]);
    cmul(cfTemp2[i],cfTemp2[i],u);
    u.re = cfTemp2[i].re;
    u.im = cfTemp2[i].im;
    cfTemp2[i].re = (u.re/(cepAmag[i]*cepAmag[i]));
    cfTemp2[i].im = (u.im/(cepAmag[i]*cepAmag[i]));
    modGd[i] = cfTemp2[i].re;
    if (i > 1)
      modGd[nfft-i+2] = modGd[i];
  }
  /* Remove DC component of modGd */

  /*RemoveAverage(modGd,nfft,&Ave); */

  /* Reduce dynamic range of modGd if necessary */
    for (i = 1; i <= nfft; i++) {
      if (modGd[i] < 0) 
	sign = -1;
      else
	sign = 1;
      abscfTemp = fabs(modGd[i]);
      if (abscfTemp == 0)
	modGd[i] = 0;
      else {
	logcfTemp = log(abscfTemp);
	if ((alfaN == 0) && (alfaP == 0))
	  if (gdsign == 1)
	    modGd[i] = sign*logcfTemp;
	  else
	    modGd[i] = logcfTemp;
	else {
	  if (sign > 0.0)
	    modGd[i] = exp(alfaP*logcfTemp);
	  else 
	    modGd[i] = exp(alfaN*logcfTemp);
	  if (gdsign == 1)
	    modGd[i] = sign*modGd[i];
	}
      }
    }   
    /* Bandlimit modGd */
    if (removeLPhase) {
      for (i = 2; i < nfBy2; i++) {
	modGd[i] = modGd[i] - modGd[1];
	modGd[nfft - i +2] = modGd[i];
      }
      modGd[1] = 0;
    }

    if ((startIndex > 1) || (endIndex < nfBy2)) {
      for (i = 1; i < startIndex; i++) { 
	modGd[i] = 0.0;
	if (i > 1)
	  modGd[nfft - i+2] = modGd[i];
      }
      for (i = endIndex+1; i <= nfBy2; i++) {
	modGd[i] = 0;
	if (i > 1)
	  modGd[nfft - i+ 2] = 0;
      }
    }
  for (i = 1; i <= nfBy2;i++)
    fvect->array[i-1] = modGd[i];
  LinearTransformationOfFVector(fvect, discreteCosineTransform, winlen, nfBy2, modGdFvect);
  for (i = 1; i <= winlen; i++) 
    modGdCepstrum[i-1] = modGdFvect->array[i-1]/(float)i/(float)scaleDCT;
  return(modGdCepstrum);  
}

/****************************************************************************
*	The following function computes  the
*	modified group delay functions from ModGdCepstrum using DCT
*
*	inputs : signal - npts
*	mfft : fft stages nfft = 2**mfft
*	winlen : window length for zero spectrum and smoothing
*
*	Outputs :
*	modGd(nfft) - smoothed modified group delay function
*****************************************************************************
*
*/
float *ModGd_DCT(float *signal,int npts,int nfft,int mfft, 
	     int winlen, int smthWinSize, float *modGd, float alfaP, float alfaN, float gamma, int gdsign, int removeLPhase,
	     int startIndex, int endIndex, float scaleDCT) {
  static float		*ax, *ay, *modGdCepstrum;
  static int              flag = 0;
  int		        i;
  float		        c0,Ave;
  if (flag == 0) {
    ax = (float *) calloc(nfft+1, sizeof(float));
    ay = (float *) calloc(nfft+1, sizeof(float));
    modGdCepstrum = (float *) calloc(winlen, sizeof(float));
    flag = 1;
  }
  modGdCepstrum = (float *) ModGdCepstrum_DCT(signal, npts, nfft, 
					  mfft, winlen, smthWinSize, modGdCepstrum, alfaP, alfaN, gamma, 
					      gdsign, removeLPhase, startIndex, endIndex, scaleDCT);
  ax[1] = 0.0;
  for (i = 2; i <= winlen; i++) {
    ax[i] = modGdCepstrum[i-2]*(i-1);
    ax[nfft-i+2] = ax[i];
  }
  for (i = winlen+1; i <= nfft-(winlen-1); i++)
    ax[i] = 0.0;
  Rfft(ax,modGd,ay,mfft,nfft,-1);
  return(modGd);
}

/****************************************************************************
*	The following function computes  the
*	modified group delay functions for the given signal
*
*	inputs : signal - npts
*	mfft : fft stages nfft = 2**mfft
*	winlen : window length for zero spectrum and smoothing
*
*	Outputs :
*	modGd(nfft) - smoothed modified group delay function
*****************************************************************************
*
*/
float *ModGd(float *signal,int npts,int nfft,int mfft, 
	     int winlen, int smthWinSize, float *modGd, 
	     float alfaP, float alfaN, float gamma, int gdsign, 
             int removeLPhase, int removeMin, 
	     int startIndex, int endIndex) {
  static float		*ax, *ay, *modGdCepstrum;
  static int              flag = 0;
  int		        i;
  float		        c0,Ave;
  if (flag == 0) {
    ax = (float *) calloc(nfft+1, sizeof(float));
    ay = (float *) calloc(nfft+1, sizeof(float));
    modGdCepstrum = (float *) calloc(winlen, sizeof(float));
    flag = 1;
  }
  modGdCepstrum = (float *) ModGdCepstrum(signal, npts, nfft, 
					  mfft, winlen, smthWinSize,
					  modGdCepstrum, alfaP, alfaN, gamma, 
					  gdsign, removeLPhase, removeMin, 
                                          startIndex, endIndex);
  ax[1] = modGdCepstrum[0];
  for (i = 2; i <= winlen; i++) {
    ax[i] = modGdCepstrum[i-1]*(i-1);
    ax[nfft-i+2] = ax[i];
  }
  for (i = winlen+1; i <= nfft-(winlen-1); i++)
    ax[i] = 0.0;
  Rfft(ax,modGd,ay,mfft,nfft,-1);
  return(modGd);
}


/****************************************************************************
*	The following function computes  the standard
*	modified group delay functions for the given signal
*
*	inputs : signal - npts
*	mfft : fft stages nfft = 2**mfft
*	winlen : window length for zero spectrum and smoothing
*
*	Outputs :
*	modGd(nfft) - smoothed modified group delay function
*****************************************************************************
*
*/
float *SmoothModGd(float *signal,int npts,int nfft,int mfft, int winlen, int smthWinSize, float *modGd) {
  static float		      *ax, *ay, *modGdTemp, *modGdCepstrum;
  static int                  flag = 0;
  /*static VECTOR_OF_F_VECTORS  *discreteInvCosineTransform;*/
  static F_VECTOR             *fvect, *modGdFvect;
  static int                  nfBy2;
  int		              i;

   if (flag == 0) {
     nfBy2 = nfft/2;
    ax = (float *) calloc(nfft+1, sizeof(float));
    ay = (float *) calloc(nfft+1, sizeof(float));
    modGdCepstrum = (float *) calloc(winlen, sizeof(float));
    modGdTemp = (float *) calloc(nfft+1, sizeof(float));
    flag = 1;
  }
modGdCepstrum = (float *) ModGdCepstrumNcN(signal, npts, nfft, 
					  mfft, winlen, smthWinSize, modGdCepstrum, 
					  1.0, 1.0, 1.0, 1, 0, 0, 1, nfft);
  ax[1] = modGdCepstrum[0];
  for (i = 2; i <= winlen; i++) {
    ax[i] = modGdCepstrum[i-1];
    ax[nfft-i+2] = ax[i];
  }
  for (i = winlen+1; i <= nfft-(winlen-1); i++)
    ax[i] = 0.0;
  Rfft(ax,modGdTemp,ay,mfft,nfft,-1);
  for (i = 0; i < nfBy2; i++)
    modGd[i] = modGdTemp[i+1];
  return(modGd);
}


/****************************************************************************
*	The following function computes  the standard
*	modified group delay functions for the given signal
*
*	inputs : signal - npts
*	mfft : fft stages nfft = 2**mfft
*	winlen : window length for zero spectrum and smoothing
*
*	Outputs :
*	modGd(nfft) - smoothed modified group delay function
*****************************************************************************
*
*/
float *StandardModGd(float *signal,int npts,int nfft,int mfft, int smthWinSize, float *modGd) {
  static float		*sigTemp;
  static float		*ax, *ay, *amag, *phase;
  static float		*cepAmag;
  static int            nfBy2, flag = 0;
  static complex        *cSig, *cfSig;
  static complex	*cfTemp1, *cfTemp2, u;
  int		        i, sign;
  float		        c0, Ave, rmax, abscfTemp, logcfTemp, alfaLogTemp;
  
  if (flag == 0) {
    nfBy2 = nfft/2;
    ax = (float *) calloc(nfft+1, sizeof(float));
    ay = (float *) calloc(nfft+1, sizeof(float));
    amag = (float *) calloc(nfft+1, sizeof(float));
    cepAmag = (float *) calloc(nfft+1, sizeof(float));
    phase = (float *) calloc(nfft+1, sizeof(float));
    cSig = (complex *) calloc (nfft+1, sizeof(complex));
    cfSig = (complex *) calloc (nfft+1, sizeof(complex));
    cfTemp1 = (complex *) calloc (nfft+1, sizeof(complex));
    cfTemp2 = (complex *) calloc (nfft+1, sizeof(complex));
    sigTemp = (float *) calloc(nfft+1, sizeof(float));
    flag = 1;
  }
  for (i = 1; i <= npts; i++)
    sigTemp[i] = signal[i];
  for(i = npts+1; i<= nfft; i++)
    sigTemp[i] = 0.0;
  for (i = 1; i < npts; i++){
    u.re = sigTemp[i];
    u.im = 0.0;
    cSig[i] = u;
    u.re = (float)(i-1);
    u.im = 0.0;
    cmul(cfSig[i],u,cSig[i]);
  }
  for (i = npts+1; i <= nfft; i++) {
    cSig[i].re = 0.0;
    cSig[i].im = 0.0;
    cfSig[i].re = 0.0;
    cfSig[i].im = 0.0;
  }
  Cfft(cSig,cfTemp1,mfft,nfft,-1);
  Cfft(cfSig,cfTemp2,mfft,nfft,-1);
  Rfft(sigTemp,ax,ay,mfft,nfft,-1);
  SpectrumReal(nfft,ax,ay,amag,phase);
  for( i = 1; i <= nfft; i++)
    ax[i] = amag[i];
  CepSmooth(ax,cepAmag,mfft,nfft, smthWinSize, &c0, 1.0);
  for (i = 1; i <= nfBy2; i++){
    conjg(u,cfTemp1[i]);
    cmul(cfTemp2[i],cfTemp2[i],u);
    u.re = cfTemp2[i].re;
    u.im = cfTemp2[i].im;
    cfTemp2[i].re = (u.re/(cepAmag[i]*cepAmag[i]));
    cfTemp2[i].im = (u.im/(cepAmag[i]*cepAmag[i]));
    modGd[i-1] = cfTemp2[i].re;
  }
return(modGd);
}


/****************************************************************************
*	The following function computes  the ModGdCepstrum from 
*       the group delay function
*
*	inputs : signal - npts
*	mfft : fft stages nfft = 2**mfft
*	winlen : window length for zero spectrum and smoothing
*
*	Outputs :
*	modGd(nfft) - ModGd Cepstra
*****************************************************************************
*
*/
float *ModGdCepstrumFromGD(float *modGd,int nfft,
	     int winlen, float *modGdCepstrum, float alfaP, float alfaN,  
		     float gamma, int gdsign, int removeLPhase, int startIndex, int endIndex, float scaleDCT) {
  static float		     *ax, *ay;
  static int                 flag = 0, nfBy2;
  int		             i;
  int                        sign;
  float                      abscfTemp, logcfTemp;
  static VECTOR_OF_F_VECTORS *discreteCosineTransform;
  static F_VECTOR            *fvect, *modGdFvect;

  if (flag == 0) {
    nfBy2 = nfft/2;
    ax = (float *) calloc(nfft+1, sizeof(float));
    ay = (float *) calloc(nfft+1, sizeof(float));
    discreteCosineTransform = (VECTOR_OF_F_VECTORS *) GeneratePseudoDct(1, winlen, nfBy2);
    fvect = (F_VECTOR *) AllocFVector(nfBy2);
    modGdFvect = (F_VECTOR *) AllocFVector(winlen);
    flag = 1;
  }

  /* RemoveAverage(modGd,nfft,&Ave); */


  /* Reduce dynamic range of modGd if necessary */
    for (i = 0; i < nfBy2; i++) {
      if (modGd[i] < 0) 
	sign = -1;
      else
	sign = 1;
      abscfTemp = fabs(modGd[i]);
      if (abscfTemp == 0)
	modGd[i] = 0;
      else {
	logcfTemp = log(abscfTemp);
	if ((alfaN == 0) && (alfaP == 0))
	  if (gdsign == 1)
	    modGd[i] = sign*logcfTemp;
	  else
	    modGd[i] = logcfTemp;
	else {
	  if (sign > 0.0)
	    modGd[i] = exp(alfaP*logcfTemp);
	  else 
	    modGd[i] = exp(alfaN*logcfTemp);
	  if (gdsign == 1)
	    modGd[i] = sign*modGd[i];
	}
      }
    }   
    /* Bandlimit modGd */
    if (removeLPhase) {
      for (i = 2; i < nfBy2; i++) {
	modGd[i] = modGd[i] - modGd[1];
	modGd[nfft - i +2] = modGd[i];
      }
      modGd[1] = 0;
    }

    if ((startIndex > 0) || (endIndex < nfBy2-1)) {
      for (i = 1; i < startIndex; i++)  
	modGd[i] = 0.0;
      for (i = endIndex+1; i <= nfBy2; i++) 
	modGd[i] = 0;
    }
  for (i = 0; i < nfBy2;i++)
    fvect->array[i] = modGd[i];
  LinearTransformationOfFVector(fvect, discreteCosineTransform, winlen, nfBy2, modGdFvect);
  for (i = 0; i < winlen; i++) 
    modGdCepstrum[i] = modGdFvect->array[i]/(float)(i+1)/(float)scaleDCT;
  return(modGdCepstrum);  
}


/****************************************************************************
*	The following function computes  the ModGdCepstrum from 
*       the Standard modgroupdelay cepstra
*
*	inputs : signal - winlen
*	mfft : fft stages nfft = 2**mfft
*	winlen : window length for zero spectrum and smoothing
*
*	Outputs :
*	modGd(nfft) - ModGd Cepstra
*****************************************************************************
*
*/
float *ModGdCepstrumFromStdCepstra(float *inModGdCepstrum,int nfft, int mfft, 
	     int winlen, float *modGdCepstrum, float alfaP, float alfaN,  
		     float gamma, int gdsign, int removeLPhase, int startIndex, int endIndex, float scaleDCT) {
  static float		     *ax, *ay, *modGdTemp, *modGd;
  static int                 flag = 0, nfBy2;
  int		             i;
  int                        sign;
  float                      abscfTemp, logcfTemp;
  static VECTOR_OF_F_VECTORS *discreteCosineTransform;
  static F_VECTOR            *fvect, *modGdFvect;

  if (flag == 0) {
    nfBy2 = nfft/2;
    ax = (float *) calloc(nfft+1, sizeof(float));
    ay = (float *) calloc(nfft+1, sizeof(float));
    discreteCosineTransform = (VECTOR_OF_F_VECTORS *) GeneratePseudoDct(1, winlen, nfBy2);
    modGd = (float *) calloc(nfft+1, sizeof(float));
    modGdTemp = (float *) calloc(nfft+1, sizeof(float));
    fvect = (F_VECTOR *) AllocFVector(nfBy2);
    modGdFvect = (F_VECTOR *) AllocFVector(winlen);
    
    flag = 1;
  }

  ax[1] = inModGdCepstrum[0];
  for (i = 2; i <= winlen+1; i++) {
    ax[i] = inModGdCepstrum[i-1];
    ax[nfft-i+2] = ax[i];
  }
  for (i = winlen+2; i <= nfft-(winlen); i++)
    ax[i] = 0.0;
  Rfft(ax,modGdTemp,ay,mfft,nfft,-1);
  for (i = 0; i < nfBy2; i++)
    modGd[i] = modGdTemp[i+1];

  /* RemoveAverage(modGd,nfft,&Ave); */

  /* Reduce dynamic range of modGd if necessary */
    for (i = 0; i < nfBy2; i++) {
      if (modGd[i] < 0) 
	sign = -1;
      else
	sign = 1;
      abscfTemp = fabs(modGd[i]);
      if (abscfTemp == 0)
	modGd[i] = 0;
      else {
	logcfTemp = log(abscfTemp);
	if ((alfaN == 0) && (alfaP == 0))
	  if (gdsign == 1)
	    modGd[i] = sign*logcfTemp;
	  else
	    modGd[i] = logcfTemp;
	else {
	  if (sign > 0.0)
	    modGd[i] = exp(alfaP*logcfTemp);
	  else 
	    modGd[i] = exp(alfaN*logcfTemp);
	  if (gdsign == 1)
	    modGd[i] = sign*modGd[i];
	}
      }
    }   
    /* Bandlimit modGd */
    if (removeLPhase) {
      for (i = 2; i < nfBy2; i++) {
	modGd[i] = modGd[i] - modGd[1];
	modGd[nfft - i +2] = modGd[i];
      }
      modGd[1] = 0;
    }

    if ((startIndex > 0) || (endIndex < nfBy2-1)) {
      for (i = 1; i < startIndex; i++)  
	modGd[i] = 0.0;
      for (i = endIndex+1; i <= nfBy2; i++) 
	modGd[i] = 0;
    }
  for (i = 0; i < nfBy2;i++)
    fvect->array[i] = modGd[i];
  LinearTransformationOfFVector (fvect, discreteCosineTransform, winlen, nfBy2, modGdFvect);
  for (i = 0; i < winlen; i++) 
    modGdCepstrum[i] = modGdFvect->array[i]/(float)(i+1)/(float)scaleDCT;
  return(modGdCepstrum);  
}

/**************************************************************************
 *                        End of DspLibrary.c
 **************************************************************************/







/*-------------------------------------------------------------------------
 * $Log: DspLibrary.c,v $
 * Revision 1.2  2004/12/23 08:28:52  hema
 * Fixed the bug SpectrumComplex
 *
 * Revision 1.1  2004/01/14 12:21:31  hema
 * Initial revision
 *
 * Revision 1.6  2002/10/16 08:53:24  hema
 * Added Modified Group Delay Cepstra
 * Min Group Delay Cepstra
 * Group Delay Cepstra
 * modified PseudoDct to include an Offset
 *
 * Revision 1.5  2002/04/23 07:26:05  hema
 * Renamed hyphenated files
 *
 * Revision 1.4  2001/10/25 12:26:32  hema
 * Modified indentation
 *
 * Revision 1.3  2000/04/28 01:57:54  hema
 * fixed group delay outlier problem
 *
 * Revision 1.2  2000/04/27 01:59:04  hema
 * modified group delay computation
 *
 * Revision 1.1  2000/04/23 02:36:13  hema
 * Initial revision
 * Local Variables:
 * time-stamp-active: t
 * time-stamp-line-limit: 20
 * time-stamp-start: "Last modified:[ 	]+"
 * time-stamp-format: "%3a %02d-%3b-%:y %02H:%02M:%02S by %u"
 * time-stamp-end: "$"
 * End:
 *                        End of DspLibrary.c
 -------------------------------------------------------------------------*/

