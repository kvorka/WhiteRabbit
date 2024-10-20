submodule (lege_poly) roots
  implicit none; contains
  
  pure real(kind=dbl) function lege_fn(deg, x)
    integer,        intent(in) :: deg
    real(kind=dbl), intent(in) :: x
    real(kind=dbl)             :: p1, p2, fac
    integer                    :: i
    
    p1      = one
    lege_fn = x
    
    do i = 2, deg
      fac     = 2 - one / i
      p2      = lege_fn * x * fac - p1 * ( fac - 1 )
      p1      = lege_fn
      lege_fn = p2
    end do
    
  end function lege_fn
  
  module subroutine find_roots_sub(this)
    class(T_legep), intent(inout) :: this
    integer                       :: i, n, ncnt
    real(kind=dbl)                :: xincr, x, fx, y, fy, t, ft, root
    real(kind=dbl), allocatable   :: xclose(:)
    
    !!***************************************************************!!
    !!* Close to roots array holder and holder arrays.              *!!
    !!***************************************************************!!
    allocate( this%rootsweights(3,this%nLege), xclose(this%nLege) )
    
    !!***************************************************************!!
    !!* Seeking for step in the bisection method.                   *!!
    !!***************************************************************!!
    n = this%nLege
      do
        n      = 2*n
        xincr  = one / n
        ncnt   = 0
        xclose = zero
        
        x  = zero
        fx = lege_fn(2*this%nLege, x)
        
        do i = 1, n
          y  = x + xincr
          fy = lege_fn(2*this%nLege, y)
          
          if (fx*fy < zero) then
            ncnt         = ncnt+1
            xclose(ncnt) = x
          end if
          
          x  = y
          fx = fy
        end do
        
        if (ncnt == this%nLege) exit
      end do
    
    !!***************************************************************!!
    !!* Seeking for roots employing bisection.                      *!!
    !!***************************************************************!!
    !$omp parallel do private (x, fx, y, fy, t, ft, root)
    do i = 1, this%nLege
      x  = xclose(i)
      y  = xclose(i) + xincr
      
      fx = lege_fn(2*this%nLege, x)
      fy = lege_fn(2*this%nLege, y)
      
      do
        t  = ( x + y ) / 2
        ft = lege_fn(2*this%nLege, t)
        
        if ( abs(ft) < 1.0d-15 ) then
          root = t
          exit
        end if
        
        if ( fx * ft < zero ) then
          y  = t
          fy = ft
        else
          x  = t
          fx = ft
        end if
        
        if ( abs( x - y ) < 1.0d-15 ) then
          root = ( x + y ) / 2
          exit
        end if
      end do
      
      this%rootsweights(1,i) = root
      this%rootsweights(2,i) = sqrt( 1 - root**2 )
      this%rootsweights(3,i) = pi * (1-root**2) / ( this%nLege * lege_fn(2*this%nLege-1, root) )**2
    end do
    !$omp end parallel do
    
    !!***************************************************************!!
    !!* Cleaning.                                                   *!!
    !!***************************************************************!!
    deallocate( xclose )
    
  end subroutine find_roots_sub
  
end submodule roots