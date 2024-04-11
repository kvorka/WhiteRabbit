submodule (PhysicalObject) Powers
  implicit none; contains
  
  module pure function buoyancy_power_fn(this) result(power)
    class(T_physicalObject), intent(in)  :: this
    integer                              :: ir
    real(kind=dbl)                       :: power
    real(kind=dbl),          allocatable :: power_ir(:)
    complex(kind=dbl),       allocatable :: gdrho_jm(:), rvelc_jm(:)
    
    allocate( power_ir(this%nd+1) , gdrho_jm(this%jms), rvelc_jm(this%jms) )
    
    do ir = 1, this%nd+1
      call this%er_buoy_rr_jm_sub( ir, gdrho_jm )
      call this%vr_rr_jm_sub( ir, rvelc_jm )
      
      power_ir(ir) = scalproduct_fn( this%jmax, gdrho_jm, rvelc_jm )
    end do
    
    power = this%rad_grid%intV_fn( power_ir )
    
    deallocate( power_ir, gdrho_jm, rvelc_jm )
    
  end function buoyancy_power_fn
  
  module pure function viscdissip_power_fn(this) result(power)
    class(T_physicalObject), intent(in)  :: this
    integer                              :: ir
    real(kind=dbl)                       :: power
    real(kind=dbl),          allocatable :: power_ir(:)
    
    allocate( power_ir(this%nd) )
    
      do ir = 1, this%nd
        power_ir(ir)  = tensnorm2_fn( this%jmax, this%sol%deviatoric_stress_jml2_fn(ir) ) / this%visc_fn(ir) / 2
      end do
      
      power = this%rad_grid%intV_fn( power_ir )
    
    deallocate( power_ir )
    
  end function viscdissip_power_fn
  
end submodule Powers