module poly_lege
use Math
implicit none; public

contains

pure subroutine zero_poly_sub(length, arr1, arr2)
  integer,         intent(in)  :: length
  complex(real64), intent(out) :: arr1(*), arr2(*)
  integer                      :: i
  
  do concurrent ( i = 1:length )
    arr1(i) = czero
    arr2(i) = czero
  end do
  
end subroutine zero_poly_sub

pure subroutine lege_setup_sub(roots, fftLege, cosx, sinx, weight)
  real(real64), intent(in)  :: roots(*), fftLege(*)
  real(real64), intent(out) :: cosx(*), sinx(*), weight(*)
  integer                   :: i
  
  do concurrent ( i = 1:step )
    cosx(i)   = roots(i)
    sinx(i)   = sqrt(1-cosx(i)**2)
    weight(i) = fftLege(i)
  end do
  
end subroutine lege_setup_sub

pure subroutine pmm_recursion_sub(m, sinx, pmm)
  integer,      intent(in)    :: m
  real(real64), intent(in)    :: sinx(*)
  real(real64), intent(inout) :: pmm(*)
  integer                     :: i
  real(real64)                :: fac
  
  select case (m)
    
    case (0)
      do concurrent ( i = 1:step )
        pmm(i) = 1._real64
      end do
    
    case default
      fac = -sqrt( ( 2*m+1 ) / ( 2._real64 * m ) )
      
      do concurrent ( i = 1:step )
        pmm(i) = fac * sinx(i) * pmm(i)
      end do
      
  end select
  
end subroutine pmm_recursion_sub

pure subroutine pmj_setup_sub(pmm, pmj2, pmj1, pmj)
  real(real64), intent(in)  :: pmm(*)
  real(real64), intent(out) :: pmj2(*), pmj1(*), pmj(*)
  integer                   :: i
  
  do concurrent ( i = 1:step )
    pmj2(i) = 0._real64
    pmj1(i) = 0._real64
    pmj(i)  = pmm(i)
  end do
  
end subroutine pmj_setup_sub

pure subroutine pmj_recursion_sub(amj, bmj, cosx, pmj2, pmj1, pmj)
  real(real64), intent(in)    :: amj, bmj, cosx(*)
  real(real64), intent(inout) :: pmj2(*), pmj1(*), pmj(*)
  integer                     :: i
  
  do concurrent ( i = 1:step )
    pmj2(i) = pmj1(i)
    pmj1(i) = pmj(i)
    pmj(i)  = amj * cosx(i) * pmj1(i) - bmj * pmj2(i)
  end do
  
end subroutine pmj_recursion_sub

pure subroutine pmj_backward_sub(nsum, pmj, cc, legesum)
  integer,         intent(in)    :: nsum
  real(real64),    intent(in)    :: pmj(*)
  complex(real64), intent(in)    :: cc(*)
  complex(real64), intent(inout) :: legesum(*)
  integer                        :: i1, i2
  
  do concurrent ( i1=1:nsum , i2=1:step )
    legesum(i2+(i1-1)*step) = legesum(i2+(i1-1)*step) + cc(i1) * pmj(i2)
  end do
  
end subroutine pmj_backward_sub

pure subroutine pmj_forward_sub(nsum, pmj, legesum, cc)
  integer,         intent(in)    :: nsum
  real(real64),    intent(in)    :: pmj(*)
  complex(real64), intent(in)    :: legesum(*)
  complex(real64), intent(inout) :: cc(*)
  integer                        :: i1, i2
  
  do concurrent ( i1=1:nsum , i2=1:step )
    cc(i1) = cc(i1) + pmj(i2) * legesum(i2+(i1-1)*step)
  end do
  
end subroutine pmj_forward_sub

pure subroutine pmj_backward_recomb_sub(m, nsum, sym, asym, sumN, sumS)
  integer,         intent(in)  :: m
  integer,         intent(in)  :: nsum
  complex(real64), intent(in)  :: sym(*), asym(*)
  complex(real64), intent(out) :: sumN(*), sumS(*)
  integer                      :: i1, i2
  
  do concurrent ( i2=1:step, i1=1:nsum )
    sumN(i1+(i2-1)*nsum) = sym(i2+(i1-1)*step) + asym(i2+(i1-1)*step)
    sumS(i1+(i2-1)*nsum) = sym(i2+(i1-1)*step) - asym(i2+(i1-1)*step)
  end do
  
  if ( m == 0 ) then
    do concurrent ( i2=1:step, i1=1:nsum )
      sumN(i1+(i2-1)*nsum)%im = 0._real64
      sumS(i1+(i2-1)*nsum)%im = 0._real64
    end do
  end if
  
end subroutine pmj_backward_recomb_sub

pure subroutine pmj_forward_recomb_sub(m, nsum, weight, sumN, sumS, sym, asym)
  integer,         intent(in)  :: m
  integer,         intent(in)  :: nsum
  real(real64),    intent(in)  :: weight(*)
  complex(real64), intent(in)  :: sumN(*), sumS(*)
  complex(real64), intent(out) :: sym(*), asym(*)
  integer                      :: i1, i2
  
  do concurrent ( i2=1:step, i1=1:nsum )
    sym(i2+(i1-1)*step)  = weight(i2) * ( sumN(i1+(i2-1)*nsum) + sumS(i1+(i2-1)*nsum) )
    asym(i2+(i1-1)*step) = weight(i2) * ( sumN(i1+(i2-1)*nsum) - sumS(i1+(i2-1)*nsum) )
  end do
  
  if ( m == 0 ) then
    do concurrent ( i2=1:step, i1=1:nsum )
      sym(i2+(i1-1)*step)%im  = 0._real64
      asym(i2+(i1-1)*step)%im = 0._real64
    end do
  end if
  
end subroutine pmj_forward_recomb_sub

end module poly_lege