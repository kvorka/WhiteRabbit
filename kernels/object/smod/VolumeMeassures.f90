submodule (PhysicalObject) VolumeMeassures
  implicit none ; contains
  
  module pure real(kind=dbl) function nuss_fn(this)
    class(T_physicalObject), intent(in) :: this
    
    nuss_fn = c2r_fn( -this%sol%flux_fn(this%nd,1,1) ) / this%r_ud / s4pi
    
  end function nuss_fn
  
  module pure real(kind=dbl) function reynolds_fn(this, choice)
    class(T_physicalObject), intent(in)           :: this
    character(len=*),        intent(in), optional :: choice
    integer                                       :: ir
    real(kind=dbl),         allocatable           :: field_vals(:)
    complex(kind=dbl),      allocatable           :: v(:)
    
    allocate( field_vals(this%nd+1), v(this%jmv) )
    
    do ir = 1, this%nd+1
      if ( ( present(choice) ) .and. ( choice == 'convective' ) ) then
        v = this%sol%conv_velocity_jml_fn(ir)
      else
        v = this%sol%velocity_jml_fn(ir)
      end if
      
      field_vals(ir) = dotproduct_fn( this%jmax, v, v )
    end do
    
    reynolds_fn = sqrt( this%rad_grid%intV_fn( field_vals ) / this%rad_grid%volume )
    
    deallocate( field_vals, v )
    
  end function reynolds_fn
  
  module pure real(kind=dbl) function volume_heating_fn(this)
    class(T_physicalObject), intent(in) :: this
    integer                             :: ir
    real(kind=dbl),         allocatable :: field_vals(:)
    complex(kind=dbl),      allocatable :: devtens(:)
    
    allocate( field_vals(this%nd), devtens(this%jmt) )
    
    do ir = 1, this%nd
      devtens        = this%sol%deviatoric_stress_jml2_fn(ir)
      field_vals(ir) = tensproduct_fn( this%jmax, devtens, devtens )
    end do
    
    volume_heating_fn = this%rad_grid%intV_fn( field_vals ) / 2
    
    deallocate(field_vals)
    
  end function volume_heating_fn

end submodule VolumeMeassures