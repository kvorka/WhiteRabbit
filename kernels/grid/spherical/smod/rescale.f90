submodule (SphericalHarmonics) rescale
  implicit none; contains
  
  module pure integer function get_maxm_fn(this, i, i2)
    class(T_lateralGrid), intent(in) :: this
    integer,              intent(in) :: i, i2
    integer                          :: m, maxm
    
    if ( allocated(this%maxm) ) then
      maxm = this%maxm(i)
      
    else
      maxm = this%jmax2
      
      do m = 0, this%jmax2
        if ( maxval(abs(this%pmm(i:i+i2-1,m))) < this%tolm ) then
          maxm = m-1
          exit
        end if
      end do
    end if
    
    get_maxm_fn = maxm
    
  end function get_maxm_fn
  
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