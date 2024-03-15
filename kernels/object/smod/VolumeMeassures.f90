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
    
    allocate( field_vals(this%nd+1) )
      
      if ( present( choice ) ) then
        
        select case (choice)
          case ('convective')
            do ir = 1, this%nd+1
              field_vals(ir) = vnorm_fn( this%jmax, this%sol%conv_velocity_jml_fn(ir) )**2
            end do
        end select
      
      else
        
        do ir = 1, this%nd+1
          field_vals(ir) = vnorm_fn( this%jmax, this%sol%velocity_jml_fn(ir) )**2
        end do
        
      end if
      
      reynolds_fn = sqrt( this%rad_grid%intV_fn( field_vals ) / this%rad_grid%volume )
    
    deallocate(field_vals)
    
  end function reynolds_fn
  
  module pure real(kind=dbl) function volume_heating_fn(this)
    class(T_physicalObject), intent(in) :: this
    integer                             :: ir
    real(kind=dbl),         allocatable :: field_vals(:)
    
    allocate( field_vals(this%nd) )
      
      do ir = 1, this%nd
        field_vals(ir) = tnorm_fn( this%jmax, this%sol%deviatoric_stress_jml2_fn(ir) )**2
      end do
      
      volume_heating_fn = this%rad_grid%intV_fn( field_vals ) / 2
    
    deallocate(field_vals)
    
  end function volume_heating_fn

end submodule VolumeMeassures