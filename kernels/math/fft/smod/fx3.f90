submodule (Fourier_transform) fx3
  implicit none
  
  real(kind=dbl), parameter :: C31 = -0.5_dbl
  real(kind=dbl), parameter :: C32 = +0.86602540378443864676_dbl
  
  contains
  
  module pure subroutine fxzm3a(m, k, l, x, t)
    integer,        intent(in)    :: m, k, l
    real(kind=dbl), intent(in)    :: t(2,0:*)
    real(kind=dbl), intent(inout) :: x(m,2,l/3,0:2,0:k-1)
    integer                       :: i, j, ij, iv
    real(kind=dbl)                :: x0re, x0im, x1re, x1im, x2re, x2im, t1re, t1im, t2re, t2im
    
    ij = 0
    
    do j = 0, k-1
      t1re = t(1,ij  )
      t1im = t(2,ij  )            
      t2re = t(1,ij+1)
      t2im = t(2,ij+1)
      
      do concurrent ( i = 1:l/3 , iv = 1:m )
        x0re =        t1re * x(iv,1,i,1,j) - t1im * x(iv,2,i,1,j)
        x0im =        t1re * x(iv,2,i,1,j) + t1im * x(iv,1,i,1,j)
        x1re = x0re - t2re * x(iv,1,i,2,j) + t2im * x(iv,2,i,2,j)
        x1im = x0im - t2re * x(iv,2,i,2,j) - t2im * x(iv,1,i,2,j)
        
        x0re = 2 * x0re          -       x1re
        x0im = 2 * x0im          -       x1im
        x2re =     x(iv,1,i,0,j) + C31 * x0re
        x2im =     x(iv,2,i,0,j) + C31 * x0im
        
        x(iv,1,i,0,j) =     x0re +       x(iv,1,i,0,j)
        x(iv,2,i,0,j) =     x0im +       x(iv,2,i,0,j)
        x(iv,1,i,2,j) =     x2re + C32 * x1im
        x(iv,2,i,2,j) =     x2im - C32 * x1re
        x(iv,1,i,1,j) = 2 * x2re -       x(iv,1,i,2,j)
        x(iv,2,i,1,j) = 2 * x2im -       x(iv,2,i,2,j)
      end do
      
      ij = ij + 2
    end do
    
  end subroutine fxzm3a
  
  module pure subroutine fxzm3b(m, l, x)
    integer,        intent(in)    :: m, l
    real(kind=dbl), intent(inout) :: x(m,2,l/3,0:2)
    integer                       :: i, iv
    real(kind=dbl)                :: x0re, x0im, x1re, x1im, x2re, x2im
    
    do concurrent ( i = 1:l/3 , iv = 1:m )
      x1re = x(iv,1,i,1) -       x(iv,1,i,2)
      x1im = x(iv,2,i,1) -       x(iv,2,i,2)
      x0re = x(iv,1,i,1) +       x(iv,1,i,2)
      x0im = x(iv,2,i,1) +       x(iv,2,i,2)
      x2re = x(iv,1,i,0) + C31 * x0re
      x2im = x(iv,2,i,0) + C31 * x0im
      
      x(iv,1,i,0) =     x0re +       x(iv,1,i,0)
      x(iv,2,i,0) =     x0im +       x(iv,2,i,0)
      x(iv,1,i,2) =     x2re + C32 * x1im
      x(iv,2,i,2) =     x2im - C32 * x1re
      x(iv,1,i,1) = 2 * x2re -       x(iv,1,i,2)
      x(iv,2,i,1) = 2 * x2im -       x(iv,2,i,2)
    end do
    
  end subroutine fxzm3b
  
end submodule fx3