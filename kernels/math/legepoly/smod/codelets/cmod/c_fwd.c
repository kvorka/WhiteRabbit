#include <stdlib.h>
#include <complex.h>

void fwd_sum( const int *n, 
              const int *step, 
              const double *pmj, 
              const double *swork, 
              double complex *cr )
{
  
  double rcc;
  double icc;
  
  for ( int i2 = 0; i2 < (*n); i2++) {
    rcc = 0.;
    icc = 0.;
    
    #pragma omp simd
    for ( int i1 = 0; i1 < (*step); i1++) {
      rcc = rcc + pmj[i1] * swork[i1+        2*(*step)*i2];
      icc = icc + pmj[i1] * swork[i1+(*step)+2*(*step)*i2];
    }
    
    cr[i2] = cr[i2] + CMPLX(rcc, icc);
  }
  
}

void fwd_shuffle( const int *n, 
                  const int *step, 
                  const double *w, 
                  const double *cosx, 
                  double *sumN, 
                  double *sumS, 
                  double *swork )
{
  double *legesumS;
  double *legesumA;

  double *sumN1;
  double *sumN2;
  double *sumS1;
  double *sumS2;
  
  for ( int i2 = 0; i2 < (*n); i2++) {
    legesumA = &swork[2*(*step)*i2               ];
    legesumS = &swork[2*(*step)*i2+2*(*step)*(*n)];
    
    sumN1 = &sumN[(*step)*i2             ];
    sumN2 = &sumN[(*step)*i2+(*n)*(*step)];
    sumS1 = &sumS[(*step)*i2             ];
    sumS2 = &sumS[(*step)*i2+(*n)*(*step)];
    
    #pragma omp simd
    for ( int i1 = 0; i1 < (*step); i1++) {
      legesumA[i1        ] = ( sumN1[i1] - sumS1[i1] ) * w[i1];
      legesumS[i1        ] = ( sumN1[i1] + sumS1[i1] ) * w[i1] * cosx[i1];
      legesumA[i1+(*step)] = ( sumN2[i1] - sumS2[i1] ) * w[i1];
      legesumS[i1+(*step)] = ( sumN2[i1] + sumS2[i1] ) * w[i1] * cosx[i1];
    }
  }
  
}