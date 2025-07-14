submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_sub
    integer :: i1, i2
    
    !$omp simd collapse (2)
    do i2 = 1, n
      do i1 = 1, step
        cr(i2) = cr(i2) + pmj(i1) * swork(i1,i2)
      end do
    end do
    
  end procedure fwd_sum_sub
  
end submodule fwd_sum
