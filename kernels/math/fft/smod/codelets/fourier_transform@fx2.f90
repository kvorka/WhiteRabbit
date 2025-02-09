submodule (fourier_transform) fx2
  implicit none; contains
  
  module procedure fxzm2a
    integer        :: i, j, iv
    real(kind=dbl) :: x1re, x1im, t1re, t1im
    
    do j = 0, k-1
      t1re = t(1,j)
      t1im = t(2,j)
      
      !$omp simd collapse(2)
      do i = 1, l/2
        do iv = 1, m
          x1re = x(iv,1,i,0,j) - t1re * x(iv,1,i,1,j)
          x1im = x(iv,2,i,0,j) - t1im * x(iv,1,i,1,j)
          
          x(iv,1,i,1,j) =     x1re          + t1im * x(iv,2,i,1,j)
          x(iv,2,i,1,j) =     x1im          - t1re * x(iv,2,i,1,j)
          x(iv,1,i,0,j) = 2 * x(iv,1,i,0,j) -        x(iv,1,i,1,j)
          x(iv,2,i,0,j) = 2 * x(iv,2,i,0,j) -        x(iv,2,i,1,j)
        end do
      end do
    end do
    
  end procedure fxzm2a
  
  module procedure fxzm2b
    integer :: i, iv
    
    !$omp simd collapse(2)
    do i = 1, l/2
      do iv = 1, m
        x(iv,1,i,1) =     x(iv,1,i,0) - x(iv,1,i,1)
        x(iv,2,i,1) =     x(iv,2,i,0) - x(iv,2,i,1)
        x(iv,1,i,0) = 2 * x(iv,1,i,0) - x(iv,1,i,1)
        x(iv,2,i,0) = 2 * x(iv,2,i,0) - x(iv,2,i,1)
      end do
    end do
    
  end procedure fxzm2b
  
end submodule fx2