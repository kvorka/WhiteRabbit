submodule (SphericalHarmonics) lege_rec_4
  implicit none; contains
  
  module pure subroutine pmj_backward_rec_4_sub(this, mj, nsum, cosx, pmj2, pmj1, pmj, cc, legesum)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: mj, nsum
    real(kind=dbl),       intent(in)    :: cosx(*)
    real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(*)
    complex(kind=dbl),    intent(inout) :: legesum(4,*)
    integer                             :: i1, i2
    
    pmj2(1:4) = pmj1(1:4)
    pmj1(1:4) = pmj(1:4)
    
    do concurrent ( i2=1:4 )
      pmj(i2)  = this%amjrr(mj) * cosx(i2) * pmj1(i2) - this%bmjrr(mj) * pmj2(i2)
    end do
    
    do concurrent ( i1=1:nsum, i2=1:4 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * pmj(i2)
    end do
    
  end subroutine pmj_backward_rec_4_sub
  
  module pure subroutine pmj_forward_rec_4_sub(this, mj, nsum, cosx, pmj2, pmj1, pmj, legesum, cr)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: mj, nsum
    real(kind=dbl),       intent(in)    :: cosx(*)
    real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: legesum(4,*)
    complex(kind=dbl),    intent(inout) :: cr(*)
    integer                             :: i1, i2
    
    pmj2(1:4) = pmj1(1:4)
    pmj1(1:4) = pmj(1:4)
    
    do concurrent ( i2=1:4 )
      pmj(i2)  = this%amjrr(mj) * cosx(i2) * pmj1(i2) - this%bmjrr(mj) * pmj2(i2)
    end do
    
    do concurrent ( i1=1:nsum , i2=1:4 )
      cr(i1) = cr(i1) + pmj(i2) * legesum(i2,i1)
    end do
    
  end subroutine pmj_forward_rec_4_sub
  
end submodule lege_rec_4