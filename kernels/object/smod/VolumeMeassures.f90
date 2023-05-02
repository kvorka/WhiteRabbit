submodule (PhysicalObject) VolumeMeassures
  implicit none
  contains
  
  real(kind=dbl) function nuss_fn(this)
    class(T_physicalObject), intent(in) :: this
    
    nuss_fn = real(-this%sol%flux_fn(this%nd,0,0,+1), kind=dbl) / this%r_ud / sqrt(4*pi)
    
  end function nuss_fn
  
  real(kind=dbl) function reynolds_fn(this)
    class(T_physicalObject), intent(in) :: this
    integer                             :: i
    real(kind=dbl),         allocatable :: field_vals(:)
    
    allocate( field_vals(this%nd+1) )
      do i = 1, this%nd+1
        field_vals(i) = vnorm_fn( this%jmax, this%sol%velocity_jml_fn(i) )**2
      end do
      
      reynolds_fn = sqrt( this%rad_grid%intV_fn( field_vals ) / this%rad_grid%volume )
    
    deallocate(field_vals)
    
  end function reynolds_fn
  
  real(kind=dbl) function nonzon_reynolds_fn(this)
    class(T_physicalObject), intent(in) :: this
    integer                             :: i
    real(kind=dbl),         allocatable :: field_vals(:)
    
    allocate( field_vals(this%nd+1) )
      do i = 1, this%nd+1
        field_vals(i) = vnorm_fn( this%jmax, this%sol%conv_velocity_jml_fn(i) )**2
      end do
      
      nonzon_reynolds_fn = sqrt( this%rad_grid%intV_fn( field_vals ) / this%rad_grid%volume )
    
    deallocate(field_vals)
    
  end function nonzon_reynolds_fn
  
  real(kind=dbl) function volume_heating_fn(this)
    class(T_physicalObject), intent(in) :: this
    integer                             :: i
    real(kind=dbl),         allocatable :: field_vals(:)
    
    allocate( field_vals(this%nd) )
      do i = 1, this%nd
        field_vals(i) = tnorm_fn( this%jmax, this%sol%deviatoric_stress_jml2_fn(i) )**2
      end do
      
      volume_heating_fn = this%rad_grid%intV_fn( field_vals ) / 2
    
    deallocate(field_vals)
    
  end function volume_heating_fn

end submodule VolumeMeassures