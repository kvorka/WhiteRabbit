module Nulify
  use Math
  implicit none; public; contains
  
  pure subroutine zero_rarray_sub( length, arr )
    integer,        intent(in)  :: length
    real(kind=dbl), intent(out) :: arr(*)
    integer                     :: i
    
    do concurrent ( i = 1:length )
      arr(i) = zero
    end do
    
  end subroutine zero_rarray_sub
  
  pure subroutine zero_carray_sub( length, arr )
    integer,           intent(in)  :: length
    complex(kind=dbl), intent(out) :: arr(*)
    integer                        :: i
    
    do concurrent ( i = 1:length )
      arr(i) = czero
    end do
    
  end subroutine zero_carray_sub
  
end module Nulify
