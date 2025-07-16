submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_sub
    integer :: i1, i2, i3
    
    do i3 = 1, n
      do i2 = 1, 2
        !$omp simd
        do i1 = 1, 16
          swork(i1,i2,i3,1) = ( sumN(i1,i3,i2) - sumS(i1,i3,i2) ) * w(i1)
          swork(i1,i2,i3,2) = ( sumN(i1,i3,i2) + sumS(i1,i3,i2) ) * w(i1) * cosx(i1)
        end do
      end do
    end do
    
  end procedure fwd_shuffle_sub
  
end submodule fwd_shuffle