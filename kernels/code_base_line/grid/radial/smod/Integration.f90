submodule(RadialGrid) Integration
  implicit none

  contains

  pure real(kind=dbl) function radial_integral_real_fn(this, field)
    class(T_radialGrid), intent(in) :: this
    real(kind=dbl),      intent(in) :: field(:)
    integer                         :: i
    real(kind=dbl), allocatable     :: field_help(:)

    allocate(field_help(this%nd))
    
      if (size(field) == this%nd+1) then
        do i = 1, this%nd
          field_help(i) = this%c(i,-1)*field(i) + this%c(i,+1)*field(i+1)
        end do
      end if

      radial_integral_real_fn = 0._dbl
        do i = 2, this%nd
          radial_integral_real_fn = radial_integral_real_fn+ (this%r(i)-this%r(i-1))*(field_help(i)+field_help(i-1))/2
        end do

    deallocate(field_help)
    
  end function radial_integral_real_fn
  
  pure real(kind=dbl) function volumetric_integral_real_fn(this, field)
    class(T_radialGrid), intent(in) :: this
    real(kind=dbl),      intent(in) :: field(:)
    integer                         :: i
    real(kind=dbl),     allocatable :: field_help(:)
    
    allocate( field_help(this%nd) )
    
      if (size(field) == this%nd+1) then
        do i = 1, this%nd
          field_help(i) = this%c(i,-1)*(this%rr(i)**2)*field(i) + this%c(i,+1)*(this%rr(i+1)**2)*field(i+1)
        end do

      else
        field_help = field * (this%r**2)
      end if

      volumetric_integral_real_fn = 0._dbl
        do i = 2, this%nd
          volumetric_integral_real_fn = volumetric_integral_real_fn + (this%r(i)-this%r(i-1))*(field_help(i)+field_help(i-1))/2
        end do

    deallocate(field_help)

  end function volumetric_integral_real_fn
  
  pure complex(kind=dbl) function radial_integral_cmplx_fn(this, field)
    class(T_radialGrid), intent(in) :: this
    complex(kind=dbl),   intent(in) :: field(:)
    integer                         :: i
    complex(kind=dbl),  allocatable :: field_help(:)

    allocate(field_help(this%nd))
    
      if (size(field) == this%nd+1) then
        do i = 1, this%nd
          field_help(i) = this%c(i,-1)*field(i) + this%c(i,+1)*field(i+1)
        end do
      end if

      radial_integral_cmplx_fn = cmplx(0._dbl, 0._dbl, kind=dbl)
        do i = 2, this%nd
          radial_integral_cmplx_fn = radial_integral_cmplx_fn + (this%r(i)-this%r(i-1))*(field_help(i)+field_help(i-1))/2
        end do

    deallocate(field_help)
    
  end function radial_integral_cmplx_fn

  pure complex(kind=dbl) function volumetric_integral_cmplx_fn(this, field)
    class(T_radialGrid), intent(in) :: this
    complex(kind=dbl),   intent(in) :: field(:)
    integer                         :: i
    complex(kind=dbl),  allocatable :: field_help(:)

    allocate( field_help(this%nd) )
    
      if (size(field) == this%nd+1) then
        do i = 1, this%nd
          field_help(i) = this%c(i,-1)*(this%rr(i)**2)*field(i) + this%c(i,+1)*(this%rr(i+1)**2)*field(i+1)
        end do

      else
        field_help = field*(this%r**2)
      end if

      volumetric_integral_cmplx_fn = cmplx(0._dbl, 0._dbl, kind=dbl)
      do i = 2, this%nd
        volumetric_integral_cmplx_fn = volumetric_integral_cmplx_fn + (this%r(i)-this%r(i-1))*(field_help(i)+field_help(i-1))/2
      end do

    deallocate(field_help)

  end function volumetric_integral_cmplx_fn

end submodule Integration