submodule (SphericalHarmonics) lege_prep_4
  implicit none; contains
  
  module pure subroutine pmj_prep_4_sub(this, i, n, cosx, weight, sumN, sumS)
    class(T_lateralGrid), intent(in)  :: this
    integer,              intent(in)  :: i, n
    real(kind=dbl),       intent(out) :: cosx(*), weight(*)
    complex(kind=dbl),    intent(out) :: sumN(*), sumS(*)
    integer                           :: sumlength
    
    cosx(1:4)   = this%roots(i:i+3)
    weight(1:4) = this%fftLege(i:i+3)
    
    sumlength = 4*n*this%jmax3
    
    call zero_carray_sub( sumlength, sumN(1) )
    call zero_carray_sub( sumlength, sumS(1) )
    
  end subroutine pmj_prep_4_sub
  
end submodule lege_prep_4