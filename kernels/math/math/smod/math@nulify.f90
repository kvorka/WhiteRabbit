submodule (math) nulify
  implicit none; contains
  
  module pure subroutine zero_rarray_sub(length, arr)
    integer,        intent(in)  :: length
    real(kind=dbl), intent(out) :: arr(*)
    integer                     :: i
    
    do concurrent ( i = 1:length )
      arr(i) = zero
    end do
    
  end subroutine zero_rarray_sub
  
  module pure subroutine zero_carray_sub(length, arr)
    integer,           intent(in)  :: length
    complex(kind=dbl), intent(out) :: arr(*)
    integer                        :: i
    
    do concurrent ( i = 1:length )
      arr(i) = czero
    end do
    
  end subroutine zero_carray_sub
  
end submodule nulify
