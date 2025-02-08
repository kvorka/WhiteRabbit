#include <stdlib.h>

void mmset( const int *ma, 
            const int *step, 
            const double *cff, 
            const double *cosx, 
            const double *sinx, 
            double *pmm, 
            double *pmj2, 
            double *pmj1, 
            double *pmj )
{
  
  if ( *ma == 1 ) {
    
    #pragma omp simd
    for ( int i2 = 0; i2 < (*step); i2++ ) {
      pmm[i2] = (*cff);
    }
    
  } else {
    
    #pragma omp simd
    for ( int i2 = 0; i2 < (*step); i2++ ) {
      pmm[i2] = (*cff) * sinx[i2] * pmm[i2];
    }
    
  }
  
  #pragma omp simd
  for ( int i2 = 0; i2 < (*step); i2++ ) {
    pmj2[i2] = 0.;
    pmj1[i2] = 0.;
    pmj[i2]  = pmm[i2] / cosx[i2];
  }
  
}

void mjrec( const int *step, 
            const double *cff, 
            const double *cosx2, 
            double *pmj2, 
            double *pmj1, 
            double *pmj )
{
  
  #pragma omp simd
  for ( int i2 = 0; i2 < (*step); i2++ ) {
    pmj2[i2] = cff[2] * pmj1[i2];
    pmj1[i2] = pmj[i2];
    pmj[i2]  = ( cff[0] * cosx2[i2] - cff[1] ) * pmj1[i2] - pmj2[i2];
  }
  
}