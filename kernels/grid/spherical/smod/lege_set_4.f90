submodule (SphericalHarmonics) lege_set_4
  implicit none; contains
  
  module pure subroutine pmj_backward_set_4_sub(this, i, m, nsum, pmj2, pmj1, pmj, cc, legesum)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: i, m, nsum
    real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(*)
    complex(kind=dbl),    intent(inout) :: legesum(4,*)
    integer                             :: i1, i2
    
    pmj2(1:4) = zero
    pmj1(1:4) = zero
    pmj(1:4)  = this%pmm(i:i+3,m)
    
    do concurrent ( i1 = 1:nsum, i2 = 1:4 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * pmj(i2)
    end do
    
  end subroutine pmj_backward_set_4_sub
  
  module pure subroutine pmj_forward_set_4_sub(this, i, m, nsum, pmj2, pmj1, pmj, legesum, cr)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: i, m, nsum
    real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: legesum(4,*)
    complex(kind=dbl),    intent(inout) :: cr(*)
    integer                             :: i1, i2
    
    pmj2(1:4) = zero
    pmj1(1:4) = zero
    pmj(1:4)  = this%pmm(i:i+3,m)
    
    do concurrent ( i1=1:nsum , i2=1:4 )
      cr(i1) = cr(i1) + pmj(i2) * legesum(i2,i1)
    end do
    
  end subroutine pmj_forward_set_4_sub
  
end submodule lege_set_4