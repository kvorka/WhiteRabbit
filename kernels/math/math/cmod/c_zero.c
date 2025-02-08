#include <stdlib.h>

void memzero( const int *n,
              double *arr )
{
  
  #pragma omp simd
  for ( int i2 = 0; i2 < (*n); i2++ ) {
    arr[i2] = 0.;
  }
  
}