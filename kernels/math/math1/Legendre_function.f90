module Legendre_function
  use Math
  implicit none; public; contains
  
  pure real(kind=dbl) function xnode_fn(nL, x, y, fx, fy)
    integer,        intent(in) :: nL
    real(kind=dbl), intent(in) :: x, y, fx, fy
    real(kind=dbl)             :: x1, fx1, x2, fx2, t, ft
    
    x1  = x
    fx1 = fx
    
    x2  = y
    fx2 = fy
    
    do
      t  = (x1 + x2)/2
      ft = lege_fn(2*nL, t)
      
      if ( abs(ft) < 1.0d-15 ) then
        xnode_fn = t ; exit
      end if
      
      if (fx1*ft < zero) then
        x2  = t
        fx2 = ft
      else
        x1  = t
        fx1 = ft
      end if
      
      if ( (x2 - x1)/abs(x1) < 1.0d-15 ) then
        xnode_fn = (x1 + x2)/2 ; exit
      end if
    end do
    
  end function xnode_fn
  
  pure real(kind=dbl) function lege_fn(deg, x)
    integer,        intent(in) :: deg
    real(kind=dbl), intent(in) :: x
    real(kind=dbl)             :: p1, p2
    integer                    :: i
    
    p1 = one
    lege_fn = x
    
    do i = 2, deg
      p2      = lege_fn * x - p1 + lege_fn * x - (lege_fn * x - p1) / i
      p1      = lege_fn
      lege_fn = p2
    end do
    
  end function lege_fn
  
end module Legendre_function
