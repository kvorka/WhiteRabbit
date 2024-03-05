submodule (SphericalHarmonics) lege_rec
  implicit none; contains
  
  module pure subroutine pmj_backward_rec_2_sub(this, mj, nsum, pmj2, pmj1, pmj, cc, legesum)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: mj, nsum
    real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(*)
    complex(kind=dbl),    intent(inout) :: legesum(2,*)
    integer                             :: i1, i2
    complex(kind=dbl)                   :: cctmp
    
    pmj2(1:2) = pmj1(1:2)
    pmj1(1:2) = pmj(1:2)
    
    do concurrent ( i = 1:2 )
      pmj(i)  = this%amjrr(mj) * cosx(i) * pmj1(i) - this%bmjrr(mj) * pmj2(i)
    end do
    
    do i1 = 1, nsum
      cctmp = cc(i1)
      
      do concurrent ( i2 = 1:2 )
        legesum(i2,i1) = legesum(i2,i1) + cctmp * pmj(i2)
      end do
    end do
    
  end subroutine pmj_backward_rec_2_sub
  
end submodule lege_rec