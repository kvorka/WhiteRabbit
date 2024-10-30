submodule (radial_grid) cvolint
  implicit none; contains
  
  module procedure volumetric_integral_cmplx_fn
    integer                        :: i
    complex(kind=dbl), allocatable :: field_help(:)
    
    allocate( field_help(2:this%nd) )
      
      if ( size(field) == this%nd+1 ) then
        !$omp parallel do
        do i = 2, this%nd
          field_help(i) = (( this%c(i-1,-1) * this%rr(i-1)**2 * field(i-1) + this%c(i-1,+1) * this%rr(i  )**2 * field(i  ) ) + &
                          &( this%c(i  ,-1) * this%rr(i  )**2 * field(i  ) + this%c(i  ,+1) * this%rr(i+1)**2 * field(i+1) )   ) / 2
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do i = 2, this%nd
          field_help(i) = ( field(i-1) * this%r(i-1)**2 + field(i) * this%r(i)**2 ) / 2
        end do
        !$omp end parallel do
      end if
      
      volumetric_integral_cmplx_fn = czero
        !$omp parallel do reduction (+:volumetric_integral_cmplx_fn)
        do i = 2, this%nd
          volumetric_integral_cmplx_fn = volumetric_integral_cmplx_fn + ( this%r(i) - this%r(i-1) ) * field_help(i)
        end do
        !$omp end parallel do
    
    deallocate( field_help )
    
  end procedure volumetric_integral_cmplx_fn
  
end submodule cvolint