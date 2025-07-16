submodule (fourier_transform) fx2
  implicit none; contains
  
  module procedure fxzm2a
    integer        :: i1, i2, i3, i4
    real(kind=dbl) :: x1re, x1im, t1re, t1im
    
    do i4 = 0, k-1
      t1re = t(1,i4)
      t1im = t(2,i4)
      
      do i3 = 1, l/2
        do i2 = 1, m
          !$omp simd
          do i1 = 1, 16
            x1re = x(i1,i2,1,i3,0,i4) - t1re * x(i1,i2,1,i3,1,i4)
            x1im = x(i1,i2,2,i3,0,i4) - t1im * x(i1,i2,1,i3,1,i4)
            
            x(i1,i2,1,i3,1,i4) =     x1re               + t1im * x(i1,i2,2,i3,1,i4)
            x(i1,i2,2,i3,1,i4) =     x1im               - t1re * x(i1,i2,2,i3,1,i4)
            x(i1,i2,1,i3,0,i4) = 2 * x(i1,i2,1,i3,0,i4) -        x(i1,i2,1,i3,1,i4)
            x(i1,i2,2,i3,0,i4) = 2 * x(i1,i2,2,i3,0,i4) -        x(i1,i2,2,i3,1,i4)
          end do
        end do
      end do
    end do
    
  end procedure fxzm2a
  
  module procedure fxzm2b
    integer :: i1, i2, i3
    
    do i3 = 1, l/2
      do i2 = 1, m
        !$omp simd
        do i1 = 1, 16
          x(i1,i2,1,i3,1) =     x(i1,i2,1,i3,0) - x(i1,i2,1,i3,1)
          x(i1,i2,2,i3,1) =     x(i1,i2,2,i3,0) - x(i1,i2,2,i3,1)
          x(i1,i2,1,i3,0) = 2 * x(i1,i2,1,i3,0) - x(i1,i2,1,i3,1)
          x(i1,i2,2,i3,0) = 2 * x(i1,i2,2,i3,0) - x(i1,i2,2,i3,1)
        end do
      end do
    end do
    
  end procedure fxzm2b
  
end submodule fx2
