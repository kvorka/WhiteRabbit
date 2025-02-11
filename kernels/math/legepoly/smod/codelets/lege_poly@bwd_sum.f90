submodule (lege_poly) bwd_sum
  implicit none; contains
  
  module procedure bwd_sum_sub
    integer        :: i1, i2
    real(kind=dbl) :: rcc, icc
    
    do concurrent ( i2 = 1:n )
      rcc = cc(i2)%re
      icc = cc(i2)%im
      
      !$omp simd
      do i1 = 1, step
        swork(i1,1,i2) = swork(i1,1,i2) + pmj(i1) * rcc
        swork(i1,2,i2) = swork(i1,2,i2) + pmj(i1) * icc
      end do
    end do
    
  end procedure bwd_sum_sub

end submodule bwd_sum
