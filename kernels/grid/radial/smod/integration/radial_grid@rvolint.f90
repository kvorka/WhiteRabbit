submodule (radial_grid) rvolint
  implicit none; contains
  
  module procedure volumetric_integral_real_fn
    integer                     :: i
    real(kind=dbl), allocatable :: field_help(:)
    
    allocate( field_help(2:this%nd) )
      
      if ( size(field) == this%nd+1 ) then
        do i = 2, this%nd
          field_help(i) = (( this%c(i-1,-1) * this%rr(i-1)**2 * field(i-1) + this%c(i-1,+1) * this%rr(i  )**2 * field(i  ) ) + &
                          &( this%c(i  ,-1) * this%rr(i  )**2 * field(i  ) + this%c(i  ,+1) * this%rr(i+1)**2 * field(i+1) )   ) / 2
        end do
      else
        do i = 2, this%nd
          field_help(i) = ( field(i-1) * this%r(i-1)**2 + field(i) * this%r(i)**2 ) / 2
        end do
      end if
    
      volumetric_integral_real_fn = zero
        do i = 2, this%nd
          volumetric_integral_real_fn = volumetric_integral_real_fn + ( this%r(i) - this%r(i-1) ) * field_help(i)
        end do
    
    deallocate( field_help )
    
  end procedure volumetric_integral_real_fn
  
end submodule rvolint