#include "malloc.h"
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
