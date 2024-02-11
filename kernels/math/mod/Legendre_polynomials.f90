module Legendre_polynomials
  use Math
  implicit none; public; contains
  
  pure real(kind=dbl) function xnode_fn(nL, x, y, fx, fy)
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
  
  pure real(kind=dbl) function lege_fn(deg, x)
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
  
  pure subroutine zero_poly_sub(length, arr1, arr2)
    integer,           intent(in)  :: length
    complex(kind=dbl), intent(out) :: arr1(*), arr2(*)
    integer                        :: i
    
    do concurrent ( i = 1:length )
      arr1(i) = czero
      arr2(i) = czero
    end do
    
  end subroutine zero_poly_sub
  
  pure subroutine lege_setup_sub(roots, fftLege, cosx, sinx, weight)
    real(kind=dbl), intent(in)  :: roots(*), fftLege(*)
    real(kind=dbl), intent(out) :: cosx(*), sinx(*), weight(*)
    integer                     :: i
    
    do concurrent ( i = 1:step )
      cosx(i)   = roots(i)
      sinx(i)   = sqrt(1-cosx(i)**2)
      weight(i) = fftLege(i)
    end do
    
  end subroutine lege_setup_sub
  
  pure subroutine pmm_recursion_sub(m, sinx, pmm)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: sinx(*)
    real(kind=dbl), intent(inout) :: pmm(*)
    integer                       :: i
    real(kind=dbl)                :: fac
    
    select case (m)
      
      case (0)
        do concurrent ( i = 1:step )
          pmm(i) = 1._dbl
        end do
      
      case default
        fac = -sqrt( ( 2*m+1 ) / ( 2._dbl * m ) )
        
        do concurrent ( i = 1:step )
          pmm(i) = fac * sinx(i) * pmm(i)
        end do
        
    end select
    
  end subroutine pmm_recursion_sub
  
  pure subroutine pmj_setup_sub(pmm, pmj2, pmj1, pmj)
    real(kind=dbl), intent(in)  :: pmm(*)
    real(kind=dbl), intent(out) :: pmj2(*), pmj1(*), pmj(*)
    integer                     :: i
    
    do concurrent ( i = 1:step )
      pmj2(i) = 0._dbl
      pmj1(i) = 0._dbl
      pmj(i)  = pmm(i)
    end do
    
  end subroutine pmj_setup_sub
  
  pure subroutine pmj_recursion_sub(amj, bmj, cosx, pmj2, pmj1, pmj)
    real(kind=dbl), intent(in)    :: amj, bmj, cosx(*)
    real(kind=dbl), intent(inout) :: pmj2(*), pmj1(*), pmj(*)
    integer                       :: i
    
    do concurrent ( i = 1:step )
      pmj2(i) = pmj1(i)
      pmj1(i) = pmj(i)
      pmj(i)  = amj * cosx(i) * pmj1(i) - bmj * pmj2(i)
    end do
    
  end subroutine pmj_recursion_sub
  
  pure subroutine pmj_backward_sub(nsum, pmj, cc, legesum)
    integer,           intent(in)    :: nsum
    real(kind=dbl),    intent(in)    :: pmj(*)
    complex(kind=dbl), intent(in)    :: cc(*)
    complex(kind=dbl), intent(inout) :: legesum(*)
    integer                          :: i1, i2
    
    do concurrent ( i1=1:nsum , i2=1:step )
      legesum(i2+(i1-1)*step) = legesum(i2+(i1-1)*step) + cc(i1) * pmj(i2)
    end do
    
  end subroutine pmj_backward_sub
  
  pure subroutine pmj_forward_sub(nsum, pmj, legesum, cc)
    integer,           intent(in)    :: nsum
    real(kind=dbl),    intent(in)    :: pmj(*)
    complex(kind=dbl), intent(in)    :: legesum(*)
    complex(kind=dbl), intent(inout) :: cc(*)
    integer                          :: i1, i2
    
    do concurrent ( i1=1:nsum , i2=1:step )
      cc(i1) = cc(i1) + pmj(i2) * legesum(i2+(i1-1)*step)
    end do
    
  end subroutine pmj_forward_sub
  
  pure subroutine pmj_backward_recomb_sub(m, nsum, sym, asym, sumN, sumS)
    integer,           intent(in)  :: m
    integer,           intent(in)  :: nsum
    complex(kind=dbl), intent(in)  :: sym(*), asym(*)
    complex(kind=dbl), intent(out) :: sumN(*), sumS(*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:step, i1=1:nsum )
      sumN(i1+(i2-1)*nsum) = sym(i2+(i1-1)*step) + asym(i2+(i1-1)*step)
      sumS(i1+(i2-1)*nsum) = sym(i2+(i1-1)*step) - asym(i2+(i1-1)*step)
    end do
    
    if ( m == 0 ) then
      do concurrent ( i2=1:step, i1=1:nsum )
        sumN(i1+(i2-1)*nsum)%im = 0._dbl
        sumS(i1+(i2-1)*nsum)%im = 0._dbl
      end do
    end if
    
  end subroutine pmj_backward_recomb_sub
  
  pure subroutine pmj_forward_recomb_sub(m, nsum, weight, sumN, sumS, sym, asym)
    integer,           intent(in)  :: m
    integer,           intent(in)  :: nsum
    real(kind=dbl),    intent(in)  :: weight(*)
    complex(kind=dbl), intent(in)  :: sumN(*), sumS(*)
    complex(kind=dbl), intent(out) :: sym(*), asym(*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:step, i1=1:nsum )
      sym(i2+(i1-1)*step)  = weight(i2) * ( sumN(i1+(i2-1)*nsum) + sumS(i1+(i2-1)*nsum) )
      asym(i2+(i1-1)*step) = weight(i2) * ( sumN(i1+(i2-1)*nsum) - sumS(i1+(i2-1)*nsum) )
    end do
    
    if ( m == 0 ) then
      do concurrent ( i2=1:step, i1=1:nsum )
        sym(i2+(i1-1)*step)%im  = 0._dbl
        asym(i2+(i1-1)*step)%im = 0._dbl
      end do
    end if
    
  end subroutine pmj_forward_recomb_sub

end module Legendre_polynomials