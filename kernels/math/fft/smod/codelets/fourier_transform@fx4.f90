submodule (fourier_transform) fx4
  implicit none; contains
  
  module procedure fxzm4a
    integer        :: i, j, ij, iv
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, t1re, t1im, t2re, t2im
    
    do concurrent ( j = 0:k-1 )
      ij = 3 * j
      
      t1re = t(1,ij  )
      t1im = t(2,ij  )
      t2re = t(1,ij+1)
      t2im = t(2,ij+1)
      
      !$omp simd collapse (2)
      do i = 1, l/4
        do iv = 1, m
          x2re = ( x(iv,1,i,0,j) - t2re * x(iv,1,i,2,j) ) + t2im * x(iv,2,i,2,j)
          x2im = ( x(iv,2,i,0,j) - t2im * x(iv,1,i,2,j) ) - t2re * x(iv,2,i,2,j)
          x3re = ( x(iv,1,i,1,j) - t2re * x(iv,1,i,3,j) ) + t2im * x(iv,2,i,3,j)
          x3im = ( x(iv,2,i,1,j) - t2im * x(iv,1,i,3,j) ) - t2re * x(iv,2,i,3,j)
          
          x0re = 2 * x(iv,1,i,0,j) - x2re
          x0im = 2 * x(iv,2,i,0,j) - x2im
          x1re = 2 * x(iv,1,i,1,j) - x3re
          x1im = 2 * x(iv,2,i,1,j) - x3im
          
          x(iv,1,i,2,j) = ( x0re - t1re * x1re ) + t1im * x1im
          x(iv,2,i,2,j) = ( x0im - t1im * x1re ) - t1re * x1im
          x(iv,1,i,0,j) = 2 * x0re - x(iv,1,i,2,j)
          x(iv,2,i,0,j) = 2 * x0im - x(Iv,2,i,2,j)
          x(iv,1,i,1,j) = ( x2re - t1re * x3im ) - t1im * x3re
          x(iv,2,i,1,j) = ( x2im + t1re * x3re ) - t1im * x3im
          x(iv,1,i,3,j) = 2 * x2re - x(iv,1,i,1,j)
          x(iv,2,i,3,j) = 2 * x2im - x(iv,2,i,1,j)
        end do
      end do
    end do
    
  end procedure fxzm4a
  
  module procedure fxzm4b
    integer        :: i, iv
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im
    
    !$omp simd collapse (2)
    do i = 1, l/4
      do iv = 1, m
        x2re = x(iv,1,i,0) - x(iv,1,i,2)
        x2im = x(iv,2,i,0) - x(iv,2,i,2)
        x0re = x(iv,1,i,0) + x(iv,1,i,2)
        x0im = x(iv,2,i,0) + x(iv,2,i,2)
        x3re = x(iv,1,i,1) - x(iv,1,i,3)
        x3im = x(iv,2,i,1) - x(iv,2,i,3)
        x1re = x(iv,1,i,1) + x(iv,1,i,3)
        x1im = x(iv,2,i,1) + x(iv,2,i,3)
        
        x(iv,1,i,2) =     x0re - x1re
        x(iv,2,i,2) =     x0im - x1im
        x(iv,1,i,0) = 2 * x0re - x(iv,1,i,2)
        x(iv,2,i,0) = 2 * x0im - x(iv,2,i,2)       
        x(iv,1,i,1) =     x2re - x3im
        x(iv,2,i,1) =     x2im + x3re
        x(iv,1,i,3) = 2 * x2re - x(iv,1,i,1)
        x(iv,2,i,3) = 2 * x2im - x(iv,2,i,1)
      end do
    end do
    
  end procedure fxzm4b
  
end submodule fx4
