submodule (lege_poly) bwd_shuffle
  implicit none; contains
  
  module procedure bwd_shuffle_sub
    integer :: i1, i2, i3
    
    !$omp simd collapse (3)
    do i3 = 1, n
      do i2 = 1, 2
        do i1 = 1, step
          swork(i1,i2,i3,2) = swork(i1,i2,i3,2) * cosx(i1)
        end do
      end do
    end do
    
    !$omp simd collapse (3)
    do i3 = 1, 2
      do i2 = 1, n
        do i1 = 1, step
          sumN(i1,i2,i3) = swork(i1,i3,i2,2) + swork(i1,i3,i2,1)
          sumS(i1,i2,i3) = swork(i1,i3,i2,2) - swork(i1,i3,i2,1)
        end do
      end do
    end do
    
  end procedure bwd_shuffle_sub

end submodule bwd_shuffle