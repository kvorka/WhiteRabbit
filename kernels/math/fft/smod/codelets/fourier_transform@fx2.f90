submodule (fourier_transform) fx2
  implicit none; contains
  
  module procedure fxzm2a
    integer        :: i, j, iv, i1
    real(kind=dbl) :: x1re, x1im, t1re, t1im
    
    do j = 0, k-1
      t1re = t(1,j)
      t1im = t(2,j)
      
      do i = 1, l/2
        do iv = 1, m
          !$omp simd
          do i1 = 1, 16
            x1re = x(i1,iv,1,i,0,j) - t1re * x(i1,iv,1,i,1,j)
            x1im = x(i1,iv,2,i,0,j) - t1im * x(i1,iv,1,i,1,j)
            
            x(i1,iv,1,i,1,j) =     x1re             + t1im * x(i1,iv,2,i,1,j)
            x(i1,iv,2,i,1,j) =     x1im             - t1re * x(i1,iv,2,i,1,j)
            x(i1,iv,1,i,0,j) = 2 * x(i1,iv,1,i,0,j) -        x(i1,iv,1,i,1,j)
            x(i1,iv,2,i,0,j) = 2 * x(i1,iv,2,i,0,j) -        x(i1,iv,2,i,1,j)
          end do
        end do
      end do
    end do
    
  end procedure fxzm2a
  
  module procedure fxzm2b
    integer :: i, iv, i1
    
    do i = 1, l/2
      do iv = 1, m
        !$omp simd
        do i1 = 1, 16
          x(i1,iv,1,i,1) =     x(i1,iv,1,i,0) - x(i1,iv,1,i,1)
          x(i1,iv,2,i,1) =     x(i1,iv,2,i,0) - x(i1,iv,2,i,1)
          x(i1,iv,1,i,0) = 2 * x(i1,iv,1,i,0) - x(i1,iv,1,i,1)
          x(i1,iv,2,i,0) = 2 * x(i1,iv,2,i,0) - x(i1,iv,2,i,1)
        end do
      end do
    end do
    
  end procedure fxzm2b
  
end submodule fx2
