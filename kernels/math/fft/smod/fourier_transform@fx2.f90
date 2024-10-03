submodule (fourier_transform) fx2
  implicit none; contains
  
  module pure subroutine fxzm2a(m, k, l, x, t)
    integer,        intent(in)    :: m, k, l
    real(kind=dbl), intent(in)    :: t(2,0:*)
    real(kind=dbl), intent(inout) :: x(m,2,l/2,0:1,0:k-1)
    integer                       :: i, j, iv
    real(kind=dbl)                :: x1re, x1im, t1re, t1im
    
    do j = 0, k-1
      t1re = t(1,j)
      t1im = t(2,j)
      
      do concurrent ( i = 1:l/2 , iv = 1:m )
        x1re = x(iv,1,i,0,j) - t1re * x(iv,1,i,1,j)
        x1im = x(iv,2,i,0,j) - t1im * x(iv,1,i,1,j)
        
        x(iv,1,i,1,j) =     x1re          + t1im * x(iv,2,i,1,j)
        x(iv,2,i,1,j) =     x1im          - t1re * x(iv,2,i,1,j)
        x(iv,1,i,0,j) = 2 * x(iv,1,i,0,j) -        x(iv,1,i,1,j)
        x(iv,2,i,0,j) = 2 * x(iv,2,i,0,j) -        x(iv,2,i,1,j)
      end do
    end do
    
  end subroutine fxzm2a
  
  module pure subroutine fxzm2b(m, l, x)
    integer,        intent(in)    :: m, l
    real(kind=dbl), intent(inout) :: x(m,2,l/2,0:1)
    integer                       :: i, iv
    
    do concurrent ( i = 1:l/2 , iv = 1:m )
      x(iv,1,i,1) =     x(iv,1,i,0) - x(iv,1,i,1)
      x(iv,2,i,1) =     x(iv,2,i,0) - x(iv,2,i,1)
      x(iv,1,i,0) = 2 * x(iv,1,i,0) - x(iv,1,i,1)
      x(iv,2,i,0) = 2 * x(iv,2,i,0) - x(iv,2,i,1)
    end do
    
  end subroutine fxzm2b
  
end submodule fx2