submodule (SphericalHarmonics) lege_prep_2
  implicit none; contains
  
  module pure subroutine pmj_prep_2_sub(this, i, n, cosx, weight, sumN, sumS)
    class(T_lateralGrid), intent(in)  :: this
    integer,              intent(in)  :: i, n
    real(kind=dbl),       intent(out) :: cosx(*), weight(*)
    complex(kind=dbl),    intent(out) :: sumN(*), sumS(*)
    integer                           :: sumlength
    
    cosx(1:2)   = this%roots(i:i+1)
    weight(1:2) = this%fftLege(i:i+1)
    
    sumlength = 2*n*this%jmax3
    
    call zero_carray_sub( sumlength, sumN(1) )
    call zero_carray_sub( sumlength, sumS(1) )
    
  end subroutine pmj_prep_2_sub
  
end submodule lege_prep_2