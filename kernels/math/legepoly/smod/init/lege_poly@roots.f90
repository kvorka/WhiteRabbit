submodule (lege_poly) roots
  implicit none; contains
  
  pure real(kind=qbl) function lege_fn(deg, x)
    integer,        intent(in) :: deg
    real(kind=qbl), intent(in) :: x
    real(kind=qbl)             :: p1, p2, fac
    integer                    :: i
    
    p1      = qone
    lege_fn = x
    
    do i = 2, deg
      fac     = 2 - qone / i
      
      p2      = fac * ( lege_fn * x - p1 ) + p1
      p1      = lege_fn
      lege_fn = p2
    end do
    
  end function lege_fn
  
  module subroutine find_roots_sub(this)
    class(T_legep), intent(inout) :: this
    integer                       :: i, n, ncnt
    real(kind=qbl)                :: xincr, x1, fx1, x2, fx2, root, froot
    real(kind=qbl), allocatable   :: xclose(:)
    
    !!**********************************************************************!!
    !!* Close to roots array holder and holder arrays.                     *!!
    !!**********************************************************************!!
    allocate( this%rootsweights(3,this%nLege), xclose(this%nLege) )
    
    !!**********************************************************************!!
    !!* Seek for efficient stepping to use within the bisection method and *!!
    !!* starting points [xclose,xclose+xincr].                             *!!
    !!**********************************************************************!!
    n = this%nLege**2 / 4
    
    do
      n     = 6 * n / 5
      xincr = qone / n
      ncnt  = 0
      
      !$omp parallel do private (fx1, fx2)
      do i = 1, n
        fx1 = lege_fn( 2*this%nLege, (i-1) * xincr )
        fx2 = lege_fn( 2*this%nLege, (i  ) * xincr )
        
        if ( fx1 * fx2 < qzero ) then
          !$omp critical
          ncnt         = ncnt+1
          xclose(ncnt) = (i-1) * xincr
          !$omp end critical
        end if
      end do
      !$omp end parallel do
      
      if ( ncnt == this%nLege ) then
        exit
      else
        do concurrent ( i = 1:ncnt )
          xclose(i) = qzero
        end do
      end if
    end do
    
    !!**********************************************************************!!
    !!* Bisection                                                          *!!
    !!**********************************************************************!!
    !$omp parallel do private (x1,fx1,x2,fx2,root,froot)
    do i = 1, this%nLege
      x1  = xclose(i)
      fx1 = lege_fn(2*this%nLege, x1)
      
      x2  = x1+xincr
      fx2 = lege_fn(2*this%nLege, x2)
      
      do
        root  = ( x1 + x2 ) / 2
        froot = lege_fn(2*this%nLege, root)
        
        if ( abs(froot) < qeps ) then
          exit
        else if ( fx1 * froot < qzero ) then
          x2  = root
          fx2 = froot
        else
          x1  = root
          fx1 = froot
        end if
      end do
      
      this%rootsweights(1,i) = root
      this%rootsweights(2,i) = sqrt( 1 - root**2 )
      this%rootsweights(3,i) = qpi * (1-root**2) / ( this%nLege * lege_fn(2*this%nLege-1, root) )**2
    end do
    !$omp end parallel do
    
    !!**********************************************************************!!
    !!* Cleaning.                                                          *!!
    !!**********************************************************************!!
    deallocate( xclose )
    
  end subroutine find_roots_sub
  
end submodule roots