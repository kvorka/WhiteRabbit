module Legendre_polynomials
  use Math
  implicit none; public; contains
  
  pure subroutine pmm_setup_2_sub(pmm)
    real(kind=dbl), intent(out) :: pmm(*)
    integer                     :: i
    
    do concurrent ( i = 1:2 )
      pmm(i) = 1._dbl
    end do
    
  end subroutine pmm_setup_2_sub
  
  pure subroutine pmm_recursion_2_sub(m, sinx, pmm)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: sinx(*)
    real(kind=dbl), intent(inout) :: pmm(*)
    integer                       :: i
    real(kind=dbl)                :: fac
    
    fac = -sqrt( ( 2*m+1 ) / ( 2._dbl * m ) )
    
    do concurrent ( i = 1:2 )
      pmm(i) = fac * sinx(i) * pmm(i)
    end do
    
  end subroutine pmm_recursion_2_sub
  
  pure subroutine pmj_setup_2_sub(pmm, pmj2, pmj1, pmj)
    real(kind=dbl), intent(in)  :: pmm(*)
    real(kind=dbl), intent(out) :: pmj2(*), pmj1(*), pmj(*)
    integer                     :: i
    
    do concurrent ( i = 1:2 )
      pmj2(i) = zero
      pmj1(i) = zero
      pmj(i)  = pmm(i)
    end do
    
  end subroutine pmj_setup_2_sub
  
  pure subroutine pmj_recursion_2_sub(amj, bmj, cosx, pmj2, pmj1, pmj)
    real(kind=dbl), intent(in)    :: amj, bmj, cosx(*)
    real(kind=dbl), intent(inout) :: pmj2(*), pmj1(*), pmj(*)
    integer                       :: i
    
    do concurrent ( i = 1:2 )
      pmj2(i) = pmj1(i)
      pmj1(i) = pmj(i)
      pmj(i)  = amj * cosx(i) * pmj1(i) - bmj * pmj2(i)
    end do
    
  end subroutine pmj_recursion_2_sub
  
  pure subroutine pmj_backward_2_sub(nsum, pmj, cc, legesum)
    integer,           intent(in)    :: nsum
    real(kind=dbl),    intent(in)    :: pmj(*)
    complex(kind=dbl), intent(in)    :: cc(*)
    complex(kind=dbl), intent(inout) :: legesum(2,*)
    integer                          :: i1, i2
    
    do concurrent ( i1=1:nsum , i2=1:2 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * pmj(i2)
    end do
    
  end subroutine pmj_backward_2_sub
  
  pure subroutine pmj_forward_2_sub(nsum, pmj, legesum, cc)
    integer,           intent(in)    :: nsum
    real(kind=dbl),    intent(in)    :: pmj(*)
    complex(kind=dbl), intent(in)    :: legesum(2,*)
    complex(kind=dbl), intent(inout) :: cc(*)
    integer                          :: i1, i2
    
    do concurrent ( i1=1:nsum , i2=1:2 )
      cc(i1) = cc(i1) + pmj(i2) * legesum(i2,i1)
    end do
    
  end subroutine pmj_forward_2_sub
  
  pure subroutine pmj_backward_recomb_2_sub(nsum, ssym, asym, sumN, sumS)
    integer,           intent(in)  :: nsum
    complex(kind=dbl), intent(in)  :: ssym(2,*), asym(2,*)
    complex(kind=dbl), intent(out) :: sumN(nsum,*), sumS(nsum,*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:2, i1=1:nsum )
      sumN(i1,i2) = ssym(i2,i1) + asym(i2,i1)
      sumS(i1,i2) = ssym(i2,i1) - asym(i2,i1)
    end do
    
  end subroutine pmj_backward_recomb_2_sub
  
  pure subroutine pmj_forward_recomb_2_sub(nsum, weight, sumN, sumS, ssym, asym)
    integer,           intent(in)  :: nsum
    real(kind=dbl),    intent(in)  :: weight(*)
    complex(kind=dbl), intent(in)  :: sumN(nsum,*), sumS(nsum,*)
    complex(kind=dbl), intent(out) :: ssym(2,*), asym(2,*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:2, i1=1:nsum )
      ssym(i2,i1) = weight(i2) * ( sumN(i1,i2) + sumS(i1,i2) )
      asym(i2,i1) = weight(i2) * ( sumN(i1,i2) - sumS(i1,i2) )
    end do
    
  end subroutine pmj_forward_recomb_2_sub
  
  pure subroutine pmm_setup_4_sub(pmm)
    real(kind=dbl), intent(out) :: pmm(*)
    integer                     :: i
    
    do concurrent ( i = 1:4 )
      pmm(i) = 1._dbl
    end do
    
  end subroutine pmm_setup_4_sub
  
  pure subroutine pmm_recursion_4_sub(m, sinx, pmm)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: sinx(*)
    real(kind=dbl), intent(inout) :: pmm(*)
    integer                       :: i
    real(kind=dbl)                :: fac
    
    fac = -sqrt( ( 2*m+1 ) / ( 2._dbl * m ) )
    
    do concurrent ( i = 1:4 )
      pmm(i) = fac * sinx(i) * pmm(i)
    end do
    
  end subroutine pmm_recursion_4_sub
  
  pure subroutine pmj_setup_4_sub(pmm, pmj2, pmj1, pmj)
    real(kind=dbl), intent(in)  :: pmm(*)
    real(kind=dbl), intent(out) :: pmj2(*), pmj1(*), pmj(*)
    integer                     :: i
    
    do concurrent ( i = 1:4 )
      pmj2(i) = zero
      pmj1(i) = zero
      pmj(i)  = pmm(i)
    end do
    
  end subroutine pmj_setup_4_sub
  
  pure subroutine pmj_recursion_4_sub(amj, bmj, cosx, pmj2, pmj1, pmj)
    real(kind=dbl), intent(in)    :: amj, bmj, cosx(*)
    real(kind=dbl), intent(inout) :: pmj2(*), pmj1(*), pmj(*)
    integer                       :: i
    
    do concurrent ( i = 1:4 )
      pmj2(i) = pmj1(i)
      pmj1(i) = pmj(i)
      pmj(i)  = amj * cosx(i) * pmj1(i) - bmj * pmj2(i)
    end do
    
  end subroutine pmj_recursion_4_sub
  
  pure subroutine pmj_backward_4_sub(nsum, pmj, cc, legesum)
    integer,           intent(in)    :: nsum
    real(kind=dbl),    intent(in)    :: pmj(*)
    complex(kind=dbl), intent(in)    :: cc(*)
    complex(kind=dbl), intent(inout) :: legesum(4,*)
    integer                          :: i1, i2
    
    do concurrent ( i1=1:nsum , i2=1:4 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * pmj(i2)
    end do
    
  end subroutine pmj_backward_4_sub
  
  pure subroutine pmj_forward_4_sub(nsum, pmj, legesum, cc)
    integer,           intent(in)    :: nsum
    real(kind=dbl),    intent(in)    :: pmj(*)
    complex(kind=dbl), intent(in)    :: legesum(4,*)
    complex(kind=dbl), intent(inout) :: cc(*)
    integer                          :: i1, i2
    
    do concurrent ( i1=1:nsum , i2=1:4 )
      cc(i1) = cc(i1) + pmj(i2) * legesum(i2,i1)
    end do
    
  end subroutine pmj_forward_4_sub
  
  pure subroutine pmj_backward_recomb_4_sub(nsum, ssym, asym, sumN, sumS)
    integer,           intent(in)  :: nsum
    complex(kind=dbl), intent(in)  :: ssym(4,*), asym(4,*)
    complex(kind=dbl), intent(out) :: sumN(nsum,*), sumS(nsum,*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:4, i1=1:nsum )
      sumN(i1,i2) = ssym(i2,i1) + asym(i2,i1)
      sumS(i1,i2) = ssym(i2,i1) - asym(i2,i1)
    end do
    
  end subroutine pmj_backward_recomb_4_sub
  
  pure subroutine pmj_forward_recomb_4_sub(nsum, weight, sumN, sumS, ssym, asym)
    integer,           intent(in)  :: nsum
    real(kind=dbl),    intent(in)  :: weight(*)
    complex(kind=dbl), intent(in)  :: sumN(nsum,*), sumS(nsum,*)
    complex(kind=dbl), intent(out) :: ssym(4,*), asym(4,*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:4, i1=1:nsum )
      ssym(i2,i1) = weight(i2) * ( sumN(i1,i2) + sumS(i1,i2) )
      asym(i2,i1) = weight(i2) * ( sumN(i1,i2) - sumS(i1,i2) )
    end do
    
  end subroutine pmj_forward_recomb_4_sub
  
  pure subroutine pmm_setup_8_sub(pmm)
    real(kind=dbl), intent(out) :: pmm(*)
    integer                     :: i
    
    do concurrent ( i = 1:8 )
      pmm(i) = 1._dbl
    end do
    
  end subroutine pmm_setup_8_sub
  
  pure subroutine pmm_recursion_8_sub(m, sinx, pmm)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: sinx(*)
    real(kind=dbl), intent(inout) :: pmm(*)
    integer                       :: i
    real(kind=dbl)                :: fac
    
    fac = -sqrt( ( 2*m+1 ) / ( 2._dbl * m ) )
    
    do concurrent ( i = 1:8 )
      pmm(i) = fac * sinx(i) * pmm(i)
    end do
    
  end subroutine pmm_recursion_8_sub
  
  pure subroutine pmj_setup_8_sub(pmm, pmj2, pmj1, pmj)
    real(kind=dbl), intent(in)  :: pmm(*)
    real(kind=dbl), intent(out) :: pmj2(*), pmj1(*), pmj(*)
    integer                     :: i
    
    do concurrent ( i = 1:8 )
      pmj2(i) = zero
      pmj1(i) = zero
      pmj(i)  = pmm(i)
    end do
    
  end subroutine pmj_setup_8_sub
  
  pure subroutine pmj_recursion_8_sub(amj, bmj, cosx, pmj2, pmj1, pmj)
    real(kind=dbl), intent(in)    :: amj, bmj, cosx(*)
    real(kind=dbl), intent(inout) :: pmj2(*), pmj1(*), pmj(*)
    integer                       :: i
    
    do concurrent ( i = 1:8 )
      pmj2(i) = pmj1(i)
      pmj1(i) = pmj(i)
      pmj(i)  = amj * cosx(i) * pmj1(i) - bmj * pmj2(i)
    end do
    
  end subroutine pmj_recursion_8_sub
  
  pure subroutine pmj_backward_8_sub(nsum, pmj, cc, legesum)
    integer,           intent(in)    :: nsum
    real(kind=dbl),    intent(in)    :: pmj(*)
    complex(kind=dbl), intent(in)    :: cc(*)
    complex(kind=dbl), intent(inout) :: legesum(8,*)
    integer                          :: i1, i2
    
    do concurrent ( i1=1:nsum , i2=1:8 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * pmj(i2)
    end do
    
  end subroutine pmj_backward_8_sub
  
  pure subroutine pmj_forward_8_sub(nsum, pmj, legesum, cc)
    integer,           intent(in)    :: nsum
    real(kind=dbl),    intent(in)    :: pmj(*)
    complex(kind=dbl), intent(in)    :: legesum(8,*)
    complex(kind=dbl), intent(inout) :: cc(*)
    integer                          :: i1, i2
    
    do concurrent ( i1=1:nsum , i2=1:8 )
      cc(i1) = cc(i1) + pmj(i2) * legesum(i2,i1)
    end do
    
  end subroutine pmj_forward_8_sub
  
  pure subroutine pmj_backward_recomb_8_sub(nsum, ssym, asym, sumN, sumS)
    integer,           intent(in)  :: nsum
    complex(kind=dbl), intent(in)  :: ssym(8,*), asym(8,*)
    complex(kind=dbl), intent(out) :: sumN(nsum,*), sumS(nsum,*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:8, i1=1:nsum )
      sumN(i1,i2) = ssym(i2,i1) + asym(i2,i1)
      sumS(i1,i2) = ssym(i2,i1) - asym(i2,i1)
    end do
    
  end subroutine pmj_backward_recomb_8_sub
  
  pure subroutine pmj_forward_recomb_8_sub(nsum, weight, sumN, sumS, ssym, asym)
    integer,           intent(in)  :: nsum
    real(kind=dbl),    intent(in)  :: weight(*)
    complex(kind=dbl), intent(in)  :: sumN(nsum,*), sumS(nsum,*)
    complex(kind=dbl), intent(out) :: ssym(8,*), asym(8,*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:8, i1=1:nsum )
      ssym(i2,i1) = weight(i2) * ( sumN(i1,i2) + sumS(i1,i2) )
      asym(i2,i1) = weight(i2) * ( sumN(i1,i2) - sumS(i1,i2) )
    end do
    
  end subroutine pmj_forward_recomb_8_sub
  
  pure subroutine pmm_setup_16_sub(pmm)
    real(kind=dbl), intent(out) :: pmm(*)
    integer                     :: i
    
    do concurrent ( i = 1:16 )
      pmm(i) = 1._dbl
    end do
    
  end subroutine pmm_setup_16_sub
  
  pure subroutine pmm_recursion_16_sub(m, sinx, pmm)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: sinx(*)
    real(kind=dbl), intent(inout) :: pmm(*)
    integer                       :: i
    real(kind=dbl)                :: fac
    
    fac = -sqrt( ( 2*m+1 ) / ( 2._dbl * m ) )
    
    do concurrent ( i = 1:16 )
      pmm(i) = fac * sinx(i) * pmm(i)
    end do
    
  end subroutine pmm_recursion_16_sub
  
  pure subroutine pmj_setup_16_sub(pmm, pmj2, pmj1, pmj)
    real(kind=dbl), intent(in)  :: pmm(*)
    real(kind=dbl), intent(out) :: pmj2(*), pmj1(*), pmj(*)
    integer                     :: i
    
    do concurrent ( i = 1:16 )
      pmj2(i) = zero
      pmj1(i) = zero
      pmj(i)  = pmm(i)
    end do
    
  end subroutine pmj_setup_16_sub
  
  pure subroutine pmj_recursion_16_sub(amj, bmj, cosx, pmj2, pmj1, pmj)
    real(kind=dbl), intent(in)    :: amj, bmj, cosx(*)
    real(kind=dbl), intent(inout) :: pmj2(*), pmj1(*), pmj(*)
    integer                       :: i
    
    do concurrent ( i = 1:16 )
      pmj2(i) = pmj1(i)
      pmj1(i) = pmj(i)
      pmj(i)  = amj * cosx(i) * pmj1(i) - bmj * pmj2(i)
    end do
    
  end subroutine pmj_recursion_16_sub
  
  pure subroutine pmj_backward_16_sub(nsum, pmj, cc, legesum)
    integer,           intent(in)    :: nsum
    real(kind=dbl),    intent(in)    :: pmj(*)
    complex(kind=dbl), intent(in)    :: cc(*)
    complex(kind=dbl), intent(inout) :: legesum(16,*)
    integer                          :: i1, i2
    
    do concurrent ( i1=1:nsum , i2=1:16 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * pmj(i2)
    end do
    
  end subroutine pmj_backward_16_sub
  
  pure subroutine pmj_forward_16_sub(nsum, pmj, legesum, cc)
    integer,           intent(in)    :: nsum
    real(kind=dbl),    intent(in)    :: pmj(*)
    complex(kind=dbl), intent(in)    :: legesum(16,*)
    complex(kind=dbl), intent(inout) :: cc(*)
    integer                          :: i1, i2
    
    do concurrent ( i1=1:nsum , i2=1:16 )
      cc(i1) = cc(i1) + pmj(i2) * legesum(i2,i1)
    end do
    
  end subroutine pmj_forward_16_sub
  
  pure subroutine pmj_backward_recomb_16_sub(nsum, ssym, asym, sumN, sumS)
    integer,           intent(in)  :: nsum
    complex(kind=dbl), intent(in)  :: ssym(16,*), asym(16,*)
    complex(kind=dbl), intent(out) :: sumN(nsum,*), sumS(nsum,*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:16, i1=1:nsum )
      sumN(i1,i2) = ssym(i2,i1) + asym(i2,i1)
      sumS(i1,i2) = ssym(i2,i1) - asym(i2,i1)
    end do
    
  end subroutine pmj_backward_recomb_16_sub
  
  pure subroutine pmj_forward_recomb_16_sub(nsum, weight, sumN, sumS, ssym, asym)
    integer,           intent(in)  :: nsum
    real(kind=dbl),    intent(in)  :: weight(*)
    complex(kind=dbl), intent(in)  :: sumN(nsum,*), sumS(nsum,*)
    complex(kind=dbl), intent(out) :: ssym(16,*), asym(16,*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:16, i1=1:nsum )
      ssym(i2,i1) = weight(i2) * ( sumN(i1,i2) + sumS(i1,i2) )
      asym(i2,i1) = weight(i2) * ( sumN(i1,i2) - sumS(i1,i2) )
    end do
    
  end subroutine pmj_forward_recomb_16_sub
  
end module Legendre_polynomials