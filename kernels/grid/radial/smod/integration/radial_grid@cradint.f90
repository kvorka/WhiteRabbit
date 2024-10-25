submodule (radial_grid) cradint
  implicit none; contains
  
  module pure complex(kind=dbl) function radial_integral_cmplx_fn(this, field)
    class(T_radialGrid), intent(in) :: this
    complex(kind=dbl),   intent(in) :: field(:)
    integer                         :: i
    complex(kind=dbl),  allocatable :: field_help(:)
    
    allocate( field_help(2:this%nd) )
      
      if (size(field) == this%nd+1) then
        do i = 2, this%nd
          field_help(i) = ( ( this%c(i-1,-1) * field(i-1) + this%c(i-1,+1) * field(i  ) ) + &
                          & ( this%c(i  ,-1) * field(i  ) + this%c(i  ,+1) * field(i+1) )   ) / 2
        end do
      else
        do i = 2, this%nd
          field_help(i) = ( field(i-1) + field(i) ) / 2
        end do
      end if
      
      radial_integral_cmplx_fn = czero
        do i = 2, this%nd
          radial_integral_cmplx_fn = radial_integral_cmplx_fn + ( this%r(i) - this%r(i-1) ) * field_help(i)
        end do
      
    deallocate( field_help )
    
  end function radial_integral_cmplx_fn
  
end submodule cradint