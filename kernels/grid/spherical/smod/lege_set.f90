submodule (SphericalHarmonics) lege_set
  implicit none; contains
  
  module pure subroutine pmj_backward_set_2_sub(this, nsum, pmm, pmj2, pmj1, pmj, cc, legesum)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: nsum
    real(kind=dbl),       intent(in)    :: pmm(*)
    real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(*)
    complex(kind=dbl),    intent(inout) :: legesum(2,*)
    integer                             :: i1, i2
    complex(kind=dbl)                   :: cctmp
    
    pmj2(1:2) = zero
    pmj1(1:2) = zero
    pmj(1:2)  = pmm(1:2)
    
    do i1 = 1, nsum
      cctmp = cc(i1)
      
      do concurrent ( i2 = 1:2 )
        legesum(i2,i1) = legesum(i2,i1) + cctmp * pmj(i2)
      end do
    end do
    
  end subroutine pmj_backward_set_2_sub
  
end submodule lege_set