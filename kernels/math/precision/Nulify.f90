module Nulify
  use Math
  implicit none; public; contains
  
  pure subroutine zero_carray_sub( length, arr )
    integer,         intent(in)  :: length
    complex(real64), intent(out) :: arr(*)
    integer                      :: i
    
    do concurrent ( i = 1:length )
      arr(i) = czero
    end do
    
  end subroutine zero_carray_sub
  
  pure subroutine zero_carray_imagpart_sub( length, arr )
    integer,         intent(in)    :: length
    complex(real64), intent(inout) :: arr(*)
    integer                        :: i
    
    do concurrent ( i = 1:length )
      arr(i)%im = zero
    end do
    
  end subroutine zero_carray_imagpart_sub
  
end module Nulify