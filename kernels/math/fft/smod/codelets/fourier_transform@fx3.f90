submodule (fourier_transform) fx3
  implicit none
  
  real(kind=dbl), parameter :: C31 = -0.5_dbl
  real(kind=dbl), parameter :: C32 = +0.86602540378443864676_dbl
  
  contains
  
  module procedure fxzm3a
    integer        :: i, j, ij, iv, i1
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im, t1re, t1im, t2re, t2im
    
    do j = 0, k-1
      ij = 2 * j
      
      t1re = t(1,ij  )
      t1im = t(2,ij  )            
      t2re = t(1,ij+1)
      t2im = t(2,ij+1)
      
      do i = 1, l/3
        do iv = 1, m
          !$omp simd
          do i1 = 1, 16
            x0re =        t1re * x(i1,iv,1,i,1,j) - t1im * x(i1,iv,2,i,1,j)
            x0im =        t1re * x(i1,iv,2,i,1,j) + t1im * x(i1,iv,1,i,1,j)
            x1re = x0re - t2re * x(i1,iv,1,i,2,j) + t2im * x(i1,iv,2,i,2,j)
            x1im = x0im - t2re * x(i1,iv,2,i,2,j) - t2im * x(i1,iv,1,i,2,j)
            
            x0re = 2 * x0re             -       x1re
            x0im = 2 * x0im             -       x1im
            x2re =     x(i1,iv,1,i,0,j) + C31 * x0re
            x2im =     x(i1,iv,2,i,0,j) + C31 * x0im
            
            x(i1,iv,1,i,0,j) =     x0re +       x(i1,iv,1,i,0,j)
            x(i1,iv,2,i,0,j) =     x0im +       x(i1,iv,2,i,0,j)
            x(i1,iv,1,i,2,j) =     x2re + C32 * x1im
            x(i1,iv,2,i,2,j) =     x2im - C32 * x1re
            x(i1,iv,1,i,1,j) = 2 * x2re -       x(i1,iv,1,i,2,j)
            x(i1,iv,2,i,1,j) = 2 * x2im -       x(i1,iv,2,i,2,j)
          end do
        end do
      end do
    end do
    
  end procedure fxzm3a
  
  module procedure fxzm3b
    integer        :: i, iv, i1
    real(kind=dbl) :: x0re, x0im, x1re, x1im, x2re, x2im
    
    do i = 1, l/3
      do iv = 1, m
        !$omp simd
        do i1 = 1, 16
          x1re = x(i1,iv,1,i,1) -       x(i1,iv,1,i,2)
          x1im = x(i1,iv,2,i,1) -       x(i1,iv,2,i,2)
          x0re = x(i1,iv,1,i,1) +       x(i1,iv,1,i,2)
          x0im = x(i1,iv,2,i,1) +       x(i1,iv,2,i,2)
          x2re = x(i1,iv,1,i,0) + C31 * x0re
          x2im = x(i1,iv,2,i,0) + C31 * x0im
          
          x(i1,iv,1,i,0) =     x0re +       x(i1,iv,1,i,0)
          x(i1,iv,2,i,0) =     x0im +       x(i1,iv,2,i,0)
          x(i1,iv,1,i,2) =     x2re + C32 * x1im
          x(i1,iv,2,i,2) =     x2im - C32 * x1re
          x(i1,iv,1,i,1) = 2 * x2re -       x(i1,iv,1,i,2)
          x(i1,iv,2,i,1) = 2 * x2im -       x(i1,iv,2,i,2)
        end do
      end do
    end do
    
  end procedure fxzm3b
  
end submodule fx3
