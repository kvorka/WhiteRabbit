module Clebsch_Legendre
  use Math
  implicit none
  
  public :: xnode_fn, lege_fn, cleb1_fn
  
  contains
  
  module pure real(kind=dbl) function xnode_fn(nL, x, y, fx, fy)
    integer,        intent(in) :: nL
    real(kind=dbl), intent(in) :: x, y, fx, fy
    real(kind=dbl)             :: x1, fx1, x2, fx2, t, ft
    
    x1 = x; fx1 = fx
    x2 = y; fx2 = fy
    
    do
      t = (x1 + x2)/2 ; ft = lege_fn(2*nL, t)
      
      if (abs(ft) < 1.0d-15) then
        xnode_fn = t ; exit
      end if
      
      if (fx1*ft < 0._dbl) then
        x2 = t; fx2 = ft
      else
        x1 = t; fx1 = ft
      end if
      
      if ((x2 - x1)/abs(x1) < 1.0d-15) then
        xnode_fn = (x1 + x2)/2 ; exit
      end if
    end do
    
  end function xnode_fn
  
  module pure real(kind=dbl) function lege_fn(deg, x)
    integer,        intent(in) :: deg
    real(kind=dbl), intent(in) :: x
    real(kind=dbl)             :: p1, p2
    integer                    :: i
    
    p1 = 1._dbl ; lege_fn = x
    
    do i = 2, deg
      p2      = lege_fn * x - p1 + lege_fn * x - (lege_fn * x - p1) / i
      p1      = lege_fn
      lege_fn = p2
    end do
    
  end function lege_fn
  
  module pure real(kind=dbl) function cleb1_fn(j1, m1, j2, m2, j, m)
    integer, intent(in) :: j1, m1, j2, m2, j, m
    
    select case (m2)
      case (-1)
        select case (j1-j)
          case (-1) ; cleb1_fn = +sqrt( (j-m-1) * (j-m  ) / ( (2*j-1._dbl) * (  j  ) * 2 ) )
          case ( 0) ; cleb1_fn = +sqrt( (j+m+1) * (j-m  ) / ( (  j+1._dbl) * (  j  ) * 2 ) )
          case (+1) ; cleb1_fn = +sqrt( (j+m+2) * (j+m+1) / ( (  j+1._dbl) * (2*j+3) * 2 ) )
        end select
        
      case (0)
        select case (j1-j)
          case (-1) ; cleb1_fn = +sqrt( (j+m  ) * (j-m  ) / ( (2*j-1._dbl) * (  j  ) ) )
          case ( 0) ; cleb1_fn = +sqrt( (  m  ) * (  m  ) / ( (  j+1._dbl) * (  j  ) ) )
          case (+1) ; cleb1_fn = -sqrt( (j+m+1) * (j-m+1) / ( (  j+1._dbl) * (2*j+3) ) )
        end select
        
      case (+1)
        select case (j1-j)
          case (-1) ; cleb1_fn = +sqrt( (j+m-1) * (j+m  ) / ( (2*j-1._dbl) * (  j  ) * 2 ) )
          case ( 0) ; cleb1_fn = -sqrt( (j+m  ) * (j-m+1) / ( (  j+1._dbl) * (  j  ) * 2 ) )
          case (+1) ; cleb1_fn = +sqrt( (j-m+1) * (j-m+2) / ( (  j+1._dbl) * (2*j+3) * 2 ) )
        end select
    end select
    
  end function cleb1_fn
  
end module Clebsch_Legendre