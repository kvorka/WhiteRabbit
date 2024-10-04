submodule (physicalobject) powers
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
  
  module pure function bottombnd_power_fn(this) result(power)
    class(T_physicalObject), intent(in) :: this
    real(kind=dbl)                      :: power
    complex(kind=dbl), allocatable      :: rvelc_jm(:)
    
    allocate( rvelc_jm(this%jms) )
    
    call this%vr_r_jm_sub( 1, rvelc_jm )
    power = this%Rad * this%gd * this%rd**2 * scalproduct_fn(this%jmax, this%bnd%t_dn, rvelc_jm)
    
    deallocate( rvelc_jm )
    
  end function bottombnd_power_fn
  
  module pure function upperbnd_power_fn(this) result(power)
    class(T_physicalObject), intent(in) :: this
    real(kind=dbl)                      :: power
    complex(kind=dbl), allocatable      :: rvelc_jm(:)
    
    allocate( rvelc_jm(this%jms) )
    
    call this%vr_r_jm_sub( this%nd, rvelc_jm )
    power = -this%Rau * this%gu * this%ru**2 * scalproduct_fn(this%jmax, this%bnd%t_up, rvelc_jm)
    
    deallocate( rvelc_jm )
    
  end function upperbnd_power_fn
  
end submodule powers