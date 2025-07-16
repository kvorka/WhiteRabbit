module alignement
  implicit none
  
#if defined (avx512)
  integer, parameter :: alig = 64      !memory alignement: AVX512
#elif defined (avx2)
  integer, parameter :: alig = 32      !memory alignement: AVX2
#else
  integer, parameter :: alig = 16      !memory alignement: AVX
#endif
  
end module alignement