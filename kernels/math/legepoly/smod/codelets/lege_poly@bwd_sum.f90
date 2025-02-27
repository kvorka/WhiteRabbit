submodule (lege_poly) bwd_sum
  implicit none; contains
  
  module procedure bwd_sum_sub
    integer :: i1, i2
    
    do i2 = 1, n
      !$omp simd
      do i1 = 1, step
        swork(i1,i2) = swork(i1,i2) + pmj(i1) * cc(i2)
      end do
    end do
    
  end procedure bwd_sum_sub

end submodule bwd_sum
