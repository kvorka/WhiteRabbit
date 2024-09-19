submodule (PhysicalObject) VolumeMeassures
  implicit none ; contains
  
  module pure real(kind=dbl) function nuss_fn(this)
    class(T_physicalObject), intent(in) :: this
    
    nuss_fn = c2r_fn( -this%q_r_fn(this%nd,1,1) ) / ( this%r_ud * s4pi )
    
  end function nuss_fn
  
  module real(kind=dbl) function reynolds_fn(this, choice)
    class(T_physicalObject), intent(in)           :: this
    character(len=*),        intent(in), optional :: choice
    integer                                       :: ir
    real(kind=dbl),          allocatable          :: field_vals(:)
    complex(kind=dbl),       allocatable          :: velocity(:)
    
    allocate( field_vals(this%nd+1), velocity(this%jmv) )
      
      if ( ( present(choice) ) .and. ( choice == 'convective' ) ) then
        
        !$omp parallel do private(velocity)
        do ir = 1, this%nd+1
          call this%sol%conv_velocity_jml_sub( ir, velocity )
          field_vals(ir) = vectnorm2_fn( this%jmax, velocity )
        end do
        !$omp end parallel do
        
      else
        
        !$omp parallel do private(velocity)
        do ir = 1, this%nd+1
          call this%v_rr_ijml_sub( ir, velocity )
          field_vals(ir) = vectnorm2_fn( this%jmax, velocity )
        end do
        !$omp end parallel do
        
      end if
      
      reynolds_fn = sqrt( this%rad_grid%intV_fn( field_vals ) / this%rad_grid%volume )
      
    deallocate( field_vals, velocity )
    
  end function reynolds_fn
  
  module pure real(kind=dbl) function volume_heating_fn(this)
    class(T_physicalObject), intent(in)  :: this
    integer                              :: ir
    real(kind=dbl),          allocatable :: field_vals(:)
    
    allocate( field_vals(this%nd) )
    
    do ir = 1, this%nd
      field_vals(ir) = tensnorm2_fn( this%jmax, this%sol%deviatoric_stress_jml2_fn(ir) )
    end do
    
    volume_heating_fn = this%rad_grid%intV_fn( field_vals ) / 2
    
    deallocate( field_vals )
    
  end function volume_heating_fn

end submodule VolumeMeassures