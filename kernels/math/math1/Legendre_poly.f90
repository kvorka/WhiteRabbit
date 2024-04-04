module Legendre_poly
  use Math
  implicit none; contains
  
  pure subroutine pmj_mmset_4_sub(m, cmm, snx, pmm, pmj2, pmj1, pmj0)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: cmm, snx(4)
    real(kind=dbl), intent(inout) :: pmm(4)
    real(kind=dbl), intent(out)   :: pmj2(4), pmj1(4), pmj0(4)
    integer                       :: i2
    
    if ( m /= 0 ) then
      do concurrent ( i2 = 1:4 )
        pmm(i2) = cmm * snx(i2) * pmm(i2)
      end do
    else
      do concurrent ( i2 = 1:4 )
        pmm(i2) = cmm
      end do
    end if
    
    do concurrent ( i2 = 1:4 )
      pmj2(i2) = zero
      pmj1(i2) = zero
      pmj0(i2) = pmm(i2)
    end do
    
  end subroutine pmj_mmset_4_sub
  
  pure subroutine pmj_recursion_4_sub(amj, bmj, csx, pmj2, pmj1, pmj0)
    real(kind=dbl), intent(in)    :: amj, bmj, csx(4)
    real(kind=dbl), intent(inout) :: pmj2(4), pmj1(4), pmj0(4)
    integer                       :: i2
    
    do concurrent ( i2=1:4 )
      pmj2(i2) = pmj1(i2)
      pmj1(i2) = pmj0(i2)
      pmj0(i2) = amj * csx(i2) * pmj1(i2) - bmj * pmj2(i2)
    end do
    
  end subroutine pmj_recursion_4_sub
  
  pure subroutine pmj_backward_sum_4_sub(nb, legep, cc, legesum)
    integer,           intent(in)    :: nb
    real(kind=dbl),    intent(in)    :: legep(4)
    complex(kind=dbl), intent(in)    :: cc(nb)
    complex(kind=dbl), intent(inout) :: legesum(4,nb)
    integer                          :: i1, i2
    
    do concurrent ( i1 = 1:nb, i2 = 1:4 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * legep(i2)
    end do
    
  end subroutine pmj_backward_sum_4_sub
  
  pure subroutine pmj_backward_rcb_4_sub(nb, sumsym, sumasym, sumN, sumS)
    integer,           intent(in)    :: nb
    complex(kind=dbl), intent(in)    :: sumsym(4,nb), sumasym(4,nb)
    complex(kind=dbl), intent(inout) :: sumN(nb,4), sumS(nb,4)
    
    sumN = transpose( sumsym + sumasym )
    sumS = transpose( sumsym - sumasym )
    
  end subroutine pmj_backward_rcb_4_sub
  
  pure subroutine pmj_forward_sum_4_sub(nf, legep, legesum, cr)
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: legep(4)
    complex(kind=dbl), intent(in)    :: legesum(4,nf)
    complex(kind=dbl), intent(inout) :: cr(nf)
    integer                          :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + sum( legep(1:4) * legesum(1:4,i1) )
    end do
    
  end subroutine pmj_forward_sum_4_sub
  
  pure subroutine pmj_forward_rcb_4_sub(nf, w, sumN, sumS, sumsym, sumasym)
    integer,           intent(in)  :: nf
    real(kind=dbl),    intent(in)  :: w(4)
    complex(kind=dbl), intent(in)  :: sumN(nf,4), sumS(nf,4)
    complex(kind=dbl), intent(out) :: sumsym(4,nf), sumasym(4,nf)
    integer                        :: i1, i2
    
    sumsym  = transpose( sumN + sumS )
    sumasym = transpose( sumN - sumS )
    
    do concurrent ( i1 = 1:nf, i2 = 1:4 )
      sumsym(i2,i1)  = w(i2) * sumsym(i2,i1)
      sumasym(i2,i1) = w(i2) * sumasym(i2,i1)
    end do
    
  end subroutine pmj_forward_rcb_4_sub
  
  pure subroutine pmj_mmset_8_sub(m, cmm, snx, pmm, pmj2, pmj1, pmj0)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: cmm, snx(8)
    real(kind=dbl), intent(inout) :: pmm(8)
    real(kind=dbl), intent(out)   :: pmj2(8), pmj1(8), pmj0(8)
    integer                       :: i2
    
    if ( m /= 0 ) then
      do concurrent ( i2 = 1:8 )
        pmm(i2) = cmm * snx(i2) * pmm(i2)
      end do
    else
      do concurrent ( i2 = 1:8 )
        pmm(i2) = cmm
      end do
    end if
    
    do concurrent ( i2 = 1:8 )
      pmj2(i2) = zero
      pmj1(i2) = zero
      pmj0(i2) = pmm(i2)
    end do
    
  end subroutine pmj_mmset_8_sub
  
  pure subroutine pmj_recursion_8_sub(amj, bmj, csx, pmj2, pmj1, pmj0)
    real(kind=dbl), intent(in)    :: amj, bmj, csx(8)
    real(kind=dbl), intent(inout) :: pmj2(8), pmj1(8), pmj0(8)
    integer                       :: i2
    
    do concurrent ( i2=1:8 )
      pmj2(i2) = pmj1(i2)
      pmj1(i2) = pmj0(i2)
      pmj0(i2) = amj * csx(i2) * pmj1(i2) - bmj * pmj2(i2)
    end do
    
  end subroutine pmj_recursion_8_sub
  
  pure subroutine pmj_backward_sum_8_sub(nb, legep, cc, legesum)
    integer,           intent(in)    :: nb
    real(kind=dbl),    intent(in)    :: legep(8)
    complex(kind=dbl), intent(in)    :: cc(nb)
    complex(kind=dbl), intent(inout) :: legesum(8,nb)
    integer                          :: i1, i2
    
    do concurrent ( i1 = 1:nb, i2 = 1:8 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * legep(i2)
    end do
    
  end subroutine pmj_backward_sum_8_sub
  
  pure subroutine pmj_backward_rcb_8_sub(nb, sumsym, sumasym, sumN, sumS)
    integer,           intent(in)    :: nb
    complex(kind=dbl), intent(in)    :: sumsym(8,nb), sumasym(8,nb)
    complex(kind=dbl), intent(inout) :: sumN(nb,8), sumS(nb,8)
    
    sumN = transpose( sumsym + sumasym )
    sumS = transpose( sumsym - sumasym )
    
  end subroutine pmj_backward_rcb_8_sub
  
  pure subroutine pmj_forward_sum_8_sub(nf, legep, legesum, cr)
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: legep(8)
    complex(kind=dbl), intent(in)    :: legesum(8,nf)
    complex(kind=dbl), intent(inout) :: cr(nf)
    integer                          :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + sum( legep(1:8) * legesum(1:8,i1) )
    end do
    
  end subroutine pmj_forward_sum_8_sub
  
  pure subroutine pmj_forward_rcb_8_sub(nf, w, sumN, sumS, sumsym, sumasym)
    integer,           intent(in)  :: nf
    real(kind=dbl),    intent(in)  :: w(8)
    complex(kind=dbl), intent(in)  :: sumN(nf,8), sumS(nf,8)
    complex(kind=dbl), intent(out) :: sumsym(8,nf), sumasym(8,nf)
    integer                        :: i1, i2
    
    sumsym  = transpose( sumN + sumS )
    sumasym = transpose( sumN - sumS )
    
    do concurrent ( i1 = 1:nf, i2 = 1:8 )
      sumsym(i2,i1)  = w(i2) * sumsym(i2,i1)
      sumasym(i2,i1) = w(i2) * sumasym(i2,i1)
    end do
    
  end subroutine pmj_forward_rcb_8_sub
  
  pure subroutine pmj_mmset_16_sub(m, cmm, snx, pmm, pmj2, pmj1, pmj0)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: cmm, snx(16)
    real(kind=dbl), intent(inout) :: pmm(16)
    real(kind=dbl), intent(out)   :: pmj2(16), pmj1(16), pmj0(16)
    integer                       :: i2
    
    if ( m /= 0 ) then
      do concurrent ( i2 = 1:16 )
        pmm(i2) = cmm * snx(i2) * pmm(i2)
      end do
    else
      do concurrent ( i2 = 1:16 )
        pmm(i2) = cmm
      end do
    end if
    
    do concurrent ( i2 = 1:16 )
      pmj2(i2) = zero
      pmj1(i2) = zero
      pmj0(i2) = pmm(i2)
    end do
    
  end subroutine pmj_mmset_16_sub
  
  pure subroutine pmj_recursion_16_sub(amj, bmj, csx, pmj2, pmj1, pmj0)
    real(kind=dbl), intent(in)    :: amj, bmj, csx(16)
    real(kind=dbl), intent(inout) :: pmj2(16), pmj1(16), pmj0(16)
    integer                       :: i2
    
    do concurrent ( i2=1:16 )
      pmj2(i2) = pmj1(i2)
      pmj1(i2) = pmj0(i2)
      pmj0(i2) = amj * csx(i2) * pmj1(i2) - bmj * pmj2(i2)
    end do
    
  end subroutine pmj_recursion_16_sub
  
  pure subroutine pmj_backward_sum_16_sub(nb, legep, cc, legesum)
    integer,           intent(in)    :: nb
    real(kind=dbl),    intent(in)    :: legep(16)
    complex(kind=dbl), intent(in)    :: cc(nb)
    complex(kind=dbl), intent(inout) :: legesum(16,nb)
    integer                          :: i1, i2
    
    do concurrent ( i1 = 1:nb, i2 = 1:16 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * legep(i2)
    end do
    
  end subroutine pmj_backward_sum_16_sub
  
  pure subroutine pmj_backward_rcb_16_sub(nb, sumsym, sumasym, sumN, sumS)
    integer,           intent(in)    :: nb
    complex(kind=dbl), intent(in)    :: sumsym(16,nb), sumasym(16,nb)
    complex(kind=dbl), intent(inout) :: sumN(nb,16), sumS(nb,16)
    
    sumN = transpose( sumsym + sumasym )
    sumS = transpose( sumsym - sumasym )
    
  end subroutine pmj_backward_rcb_16_sub
  
  pure subroutine pmj_forward_sum_16_sub(nf, legep, legesum, cr)
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: legep(16)
    complex(kind=dbl), intent(in)    :: legesum(16,nf)
    complex(kind=dbl), intent(inout) :: cr(nf)
    integer                          :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + sum( legep(1:16) * legesum(1:16,i1) )
    end do
    
  end subroutine pmj_forward_sum_16_sub
  
  pure subroutine pmj_forward_rcb_16_sub(nf, w, sumN, sumS, sumsym, sumasym)
    integer,           intent(in)  :: nf
    real(kind=dbl),    intent(in)  :: w(16)
    complex(kind=dbl), intent(in)  :: sumN(nf,16), sumS(nf,16)
    complex(kind=dbl), intent(out) :: sumsym(16,nf), sumasym(16,nf)
    integer                        :: i1, i2
    
    sumsym  = transpose( sumN + sumS )
    sumasym = transpose( sumN - sumS )
    
    do concurrent ( i1 = 1:nf, i2 = 1:16 )
      sumsym(i2,i1)  = w(i2) * sumsym(i2,i1)
      sumasym(i2,i1) = w(i2) * sumasym(i2,i1)
    end do
    
  end subroutine pmj_forward_rcb_16_sub
  
end module Legendre_poly