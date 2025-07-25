submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_sub
    integer :: i1, i2
    
    do i2 = 1, n
      !$omp simd
      do i1 = 1, 16
        cr(i2) = cr(i2) + pmj(i1) * swork(i1,i2)
      end do
    end do
    
  end procedure fwd_sum_sub
  
end submodule fwd_sum
