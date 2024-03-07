submodule (SphericalHarmonics) rescale
  implicit none; contains
  
  module pure subroutine rescale_sub( this, arr, length )
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: length
    complex(kind=dbl),    intent(inout) :: arr(*)
    integer                             :: ijm
    
    do concurrent ( ijm = 1:length )
      arr(ijm) = this%scale * arr(ijm)
    end do
    
  end subroutine rescale_sub
  
end submodule rescale