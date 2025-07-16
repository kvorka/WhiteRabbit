submodule (fourier_transform) fx4
  implicit none; contains
  
  module procedure fxzm4a
    integer        :: i1, i2, i3, i4
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, t1re, t1im, t2re, t2im
    
    do i4 = 0, k-1
      i1 = 3 * i4
      
      t1re = t(1,i1  )
      t1im = t(2,i1  )
      t2re = t(1,i1+1)
      t2im = t(2,i1+1)
      
      do i3 = 1, l/4
        do i2 = 1, m
          !$omp simd
          do i1 = 1, 16
            x2re = ( x(i1,i2,1,i3,0,i4) - t2re * x(i1,i2,1,i3,2,i4) ) + t2im * x(i1,i2,2,i3,2,i4)
            x2im = ( x(i1,i2,2,i3,0,i4) - t2im * x(i1,i2,1,i3,2,i4) ) - t2re * x(i1,i2,2,i3,2,i4)
            x3re = ( x(i1,i2,1,i3,1,i4) - t2re * x(i1,i2,1,i3,3,i4) ) + t2im * x(i1,i2,2,i3,3,i4)
            x3im = ( x(i1,i2,2,i3,1,i4) - t2im * x(i1,i2,1,i3,3,i4) ) - t2re * x(i1,i2,2,i3,3,i4)
            
            x0re = 2 * x(i1,i2,1,i3,0,i4) - x2re
            x0im = 2 * x(i1,i2,2,i3,0,i4) - x2im
            x1re = 2 * x(i1,i2,1,i3,1,i4) - x3re
            x1im = 2 * x(i1,i2,2,i3,1,i4) - x3im
            
            x(i1,i2,1,i3,2,i4) = ( x0re - t1re * x1re ) + t1im * x1im
            x(i1,i2,2,i3,2,i4) = ( x0im - t1im * x1re ) - t1re * x1im
            x(i1,i2,1,i3,0,i4) = 2 * x0re - x(i1,i2,1,i3,2,i4)
            x(i1,i2,2,i3,0,i4) = 2 * x0im - x(i1,i2,2,i3,2,i4)
            x(i1,i2,1,i3,1,i4) = ( x2re - t1re * x3im ) - t1im * x3re
            x(i1,i2,2,i3,1,i4) = ( x2im + t1re * x3re ) - t1im * x3im
            x(i1,i2,1,i3,3,i4) = 2 * x2re - x(i1,i2,1,i3,1,i4)
            x(i1,i2,2,i3,3,i4) = 2 * x2im - x(i1,i2,2,i3,1,i4)
          end do
        end do
      end do
    end do
    
  end procedure fxzm4a
  
  module procedure fxzm4b
    integer        :: i1, i2, i3
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im
    
    do i3 = 1, l/4
      do i2 = 1, m
        !$omp simd
        do i1 = 1, 16
          x2re = x(i1,i2,1,i3,0) - x(i1,i2,1,i3,2)
          x2im = x(i1,i2,2,i3,0) - x(i1,i2,2,i3,2)
          x0re = x(i1,i2,1,i3,0) + x(i1,i2,1,i3,2)
          x0im = x(i1,i2,2,i3,0) + x(i1,i2,2,i3,2)
          x3re = x(i1,i2,1,i3,1) - x(i1,i2,1,i3,3)
          x3im = x(i1,i2,2,i3,1) - x(i1,i2,2,i3,3)
          x1re = x(i1,i2,1,i3,1) + x(i1,i2,1,i3,3)
          x1im = x(i1,i2,2,i3,1) + x(i1,i2,2,i3,3)
          
          x(i1,i2,1,i3,2) =     x0re - x1re
          x(i1,i2,2,i3,2) =     x0im - x1im
          x(i1,i2,1,i3,0) = 2 * x0re - x(i1,i2,1,i3,2)
          x(i1,i2,2,i3,0) = 2 * x0im - x(i1,i2,2,i3,2)       
          x(i1,i2,1,i3,1) =     x2re - x3im
          x(i1,i2,2,i3,1) =     x2im + x3re
          x(i1,i2,1,i3,3) = 2 * x2re - x(i1,i2,1,i3,1)
          x(i1,i2,2,i3,3) = 2 * x2im - x(i1,i2,2,i3,1)
        end do
      end do
    end do
    
  end procedure fxzm4b
  
end submodule fx4
