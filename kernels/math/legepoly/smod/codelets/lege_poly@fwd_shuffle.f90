submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_sub
    integer :: i1, i2, i3
    
    !$omp simd collapse (3)
    do i3 = 1, n
      do i2 = 1, 2
        do i1 = 1, step
          swork(i1,i2,i3,1) = sumN(i1,i3,i2) - sumS(i1,i3,i2)
          swork(i1,i2,i3,2) = sumN(i1,i3,i2) + sumS(i1,i3,i2)
        end do
      end do
    end do
    
    !$omp simd collapse (3)
    do i3 = 1, n
      do i2 = 1, 2
        do i1 = 1, step
          swork(i1,i2,i3,1) = swork(i1,i2,i3,1) * w(i1)
          swork(i1,i2,i3,2) = swork(i1,i2,i3,2) * w(i1) * cosx(i1)
        end do
      end do
    end do
    
  end procedure fwd_shuffle_sub
  
end submodule fwd_shuffle