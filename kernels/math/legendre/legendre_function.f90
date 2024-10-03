module legendre_function
  use math
  implicit none; public; contains
  
  pure real(kind=dbl) function lege_fn(deg, x)
    integer,        intent(in) :: deg
    real(kind=dbl), intent(in) :: x
    real(kind=dbl)             :: p1, p2
    integer                    :: i
    
    p1 = one
    lege_fn = x
    
    do i = 2, deg
      p2      = lege_fn * x + (lege_fn * x - p1) * ( 1 - one / i )
      p1      = lege_fn
      lege_fn = p2
    end do
    
  end function lege_fn
  
  pure real(kind=dbl) function xstep_fn(deg)
    integer, intent(in) :: deg
    integer             :: n, ncnt, i
    real(kind=dbl)      :: x, y, fx, fy
    
    n = deg
    
    do
      n        = 2*n
      xstep_fn = one/n
      ncnt     = 0
      
      x  = zero
      fx = lege_fn(2*deg, x)
      
      do i = 1, n
        y  = x + xstep_fn
        fy = lege_fn(2*deg, y)
        
        if (fx*fy < zero) then
          ncnt = ncnt+1
        end if
        
        x  = y
        fx = fy
      end do
      
      if (ncnt == deg) exit
    end do
    
  end function xstep_fn
  
  pure real(kind=dbl) function xnode_fn(deg, x, y, fx, fy)
    integer,        intent(in) :: deg
    real(kind=dbl), intent(in) :: x, y, fx, fy
    real(kind=dbl)             :: x1, fx1, x2, fx2, t, ft
    
    x1  = x
    fx1 = fx
    
    x2  = y
    fx2 = fy
    
    do
      t  = (x1 + x2)/2
      ft = lege_fn(2*deg, t)
      
      if ( abs(ft) < 1.0d-15 ) then
        xnode_fn = t
        exit
      end if
      
      if (fx1*ft < zero) then
        x2  = t
        fx2 = ft
      else
        x1  = t
        fx1 = ft
      end if
      
      if ( (x2 - x1)/abs(x1) < 1.0d-15 ) then
        xnode_fn = (x1 + x2)/2
        exit
      end if
    end do
    
  end function xnode_fn
  
  pure subroutine legef_roots_sub(deg, roots)
    integer,                     intent(in)  :: deg
    real(kind=dbl), allocatable, intent(out) :: roots(:)
    integer                                  :: i
    real(kind=dbl)                           :: xincr, x, fx, y, fy
    
    allocate( roots(deg) )
    
    i = 0
    
    xincr = xstep_fn(deg)
    x     = zero
    y     = xincr
    fx    = lege_fn(2*deg, x)
    fy    = lege_fn(2*deg, y)
    
    do
      if (fx*fy < zero) then
        i = i+1
        roots(i) = xnode_fn(deg, x, y, fx, fy)
        
        if (i == deg) exit
      end if
      
      x  = y
      y  = x + xincr
      fx = fy
      fy = lege_fn(2*deg, y)
    end do
    
  end subroutine legef_roots_sub
  
  pure subroutine legef_weights_sub(deg, roots, weights)
    integer,                     intent(in)  :: deg
    real(kind=dbl),              intent(in)  :: roots(deg)
    real(kind=dbl), allocatable, intent(out) :: weights(:)
    integer                                  :: i
    
    allocate( weights(deg) )
    
    do i = 1, deg
      weights(i) = pi * (1-roots(i)**2) / ( deg * lege_fn(2*deg-1, roots(i)) )**2
    end do
    
  end subroutine legef_weights_sub
  
end module legendre_function
