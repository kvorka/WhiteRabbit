#include <stdlib.h>
#include <complex.h>

void bwd_sum( const int *n, 
              const int *step, 
              const double *pmj, 
              const double complex *cc, 
              double *swork )
{
  
  double rcc;
  double icc;
  
  double *legesum;
  
  for ( int i2 = 0; i2 < (*n); i2++) {
    rcc = creal( cc[i2] );
    icc = cimag( cc[i2] );
    
    legesum = &swork[2*(*step)*i2];
    
    #pragma omp simd
    for ( int i1 = 0; i1 < (*step); i1++) {
      legesum[i1        ] = legesum[i1        ] + pmj[i1] * rcc;
      legesum[i1+(*step)] = legesum[i1+(*step)] + pmj[i1] * icc;
    }

  }
  
}

void bwd_shuffle( const int *n, 
                  const int *step, 
                  const double *cosx, 
                  double *swork, 
                  double *sumN, 
                  double *sumS )
{
  
  double *legesumS;
  double *legesumA;

  double *sumN1;
  double *sumN2;
  double *sumS1;
  double *sumS2;
  
  for ( int i2 = 0; i2 < (*n); i2++) {
    legesumS = &swork[2*(*step)*i2+2*(*step)*(*n)];
    
    #pragma omp simd
    for ( int i1 = 0; i1 < (*step); i1++) {
      legesumS[i1        ] = legesumS[i1        ] * cosx[i1];
      legesumS[i1+(*step)] = legesumS[i1+(*step)] * cosx[i1];
    }
  }
  
  for ( int i2 = 0; i2 < (*n); i2++) {
    legesumA = &swork[2*(*step)*i2               ];
    legesumS = &swork[2*(*step)*i2+2*(*step)*(*n)];
    
    sumN1 = &sumN[(*step)*i2             ];
    sumN2 = &sumN[(*step)*i2+(*n)*(*step)];
    sumS1 = &sumS[(*step)*i2             ];
    sumS2 = &sumS[(*step)*i2+(*n)*(*step)];
    
    #pragma omp simd
    for ( int i1 = 0; i1 < (*step); i1++) {
      sumN1[i1] = legesumS[i1        ] + legesumA[i1        ];
      sumS1[i1] = legesumS[i1        ] - legesumA[i1        ];
      sumN2[i1] = legesumS[i1+(*step)] + legesumA[i1+(*step)];
      sumS2[i1] = legesumS[i1+(*step)] - legesumA[i1+(*step)];
    }
  }
  
}