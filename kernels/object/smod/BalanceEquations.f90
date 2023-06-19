submodule (PhysicalObject) BalanceEquations
  implicit none
  contains
  
  real(kind=dbl) function laws_mech_fn(this)
    class(T_physicalObject), intent(in) :: this
    
    select case( this%mechanic_bnd )
      case( 'shape' )
        laws_mech_fn = bound_power_fn(this) / ( heating_power_fn(this) - buoyancy_power_fn(this) )
      case default  
        laws_mech_fn = buoyancy_power_fn(this) / heating_power_fn(this)
    end select
    
  end function laws_mech_fn
  
  real(kind=dbl) function laws_temp_fn(this)
    class(T_physicalObject), intent(in) :: this
    
    select case( this%thermal_bnd )
      case( 'phase' )
        laws_temp_fn = advected_heat_fn(this) / bound_flux_fn(this)
      case default  
        laws_temp_fn = real(this%sol%flux_fn(this%nd,1,1), kind=dbl) / real(this%sol%flux_fn(1,1,1), kind=dbl) / (this%r_ud**2)
    end select
    
  end function laws_temp_fn
  
  real(kind=dbl) function laws_force_fn(this, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ijm
    integer                             :: ir
    complex(kind=dbl)                   :: press_topo_d, press_topo_u, press_buoy
    complex(kind=dbl),      allocatable :: drho(:)
    
    press_topo_u = this%Rau * this%ru**2 * this%gravity%g_fn(this%ru) * this%sol%t_up(ijm)
    press_topo_d = this%Rad * this%rd**2 * this%gravity%g_fn(this%rd) * this%sol%t_dn(ijm)
    
    allocate( drho(this%nd+1) )
      do ir = 1, this%nd+1
        drho(ir) = this%Ra * this%alpha_fn(ir) * this%gravity%g_fn(this%rad_grid%rr(ir)) * this%sol%temp_fn(ir,ijm)
      end do
      
      press_buoy = this%rad_grid%intV_fn( drho )
    deallocate( drho )
    
    laws_force_fn = real( press_topo_u / ( press_topo_d + press_buoy ), kind=dbl )
    
  end function laws_force_fn
  
  real(kind=dbl) function buoyancy_power_fn(this)
    class(T_physicalObject), intent(in) :: this
    integer                             :: ir
    real(kind=dbl),         allocatable :: power_i(:)
    
    allocate( power_i(this%nd+1) )
    
      do ir = 1, this%nd+1
        power_i(ir) = dotproduct_fn( this%jmax, this%buoy_rr_jml_fn(ir), this%sol%velocity_jml_fn(ir) )
      end do
      
      buoyancy_power_fn = this%rad_grid%intV_fn( power_i )
    
    deallocate( power_i )
    
  end function buoyancy_power_fn
  
  real(kind=dbl) function heating_power_fn(this)
    class(T_physicalObject),       intent(in) :: this
    integer                                   :: ir
    real(kind=dbl),               allocatable :: power_i(:)
    
    allocate( power_i(this%nd) )
    
      do ir = 1, this%nd
        power_i(ir) = tnorm_fn(this%jmax, this%sol%deviatoric_stress_jml2_fn(ir))**2 / this%visc_fn(ir) / 2
      end do
      
      heating_power_fn = this%rad_grid%intV_fn( power_i )
    
    deallocate( power_i )
    
  end function heating_power_fn
  
  real(kind=dbl) function bound_power_fn(this)
    class(T_physicalObject), intent(in) :: this
      
    associate( rd => this%rad_grid%r(1),                         &
             & ru => this%rad_grid%r(this%nd),                   &
             & gd => this%gravity%g_fn(this%rad_grid%r(1)),      &
             & gu => this%gravity%g_fn(this%rad_grid%r(this%nd)) )
      
      bound_power_fn = this%Rad * this%gravity%g_fn(rd) * rd**2 * scalproduct_fn(this%jmax, this%sol%t_dn, this%vr_jm_fn(1)) - &
                     & this%Rau * this%gravity%g_fn(ru) * ru**2 * scalproduct_fn(this%jmax, this%sol%t_up, this%vr_jm_fn(this%nd))
      
    end associate
    
  end function bound_power_fn
  
  real(kind=dbl) function bound_flux_fn(this)
    class(T_physicalObject), intent(in) :: this
    
    bound_flux_fn = ( real(-this%sol%flux_fn(this%nd,1,1), kind=dbl) * this%rad_grid%r(this%nd)**2 - &
                    & real(-this%sol%flux_fn(      1,1,1), kind=dbl) * this%rad_grid%r(      1)**2 - &
                    & this%Ds/this%Ra * this%rad_grid%intV_fn( real(this%htide(:,1),kind=dbl) ) ) / sqrt(4*pi)
    
  end function bound_flux_fn
  
  real(kind=dbl) function advected_heat_fn(this)
    class(T_physicalObject), intent(in)  :: this
    integer                              :: ir
    real(kind=dbl),          allocatable :: heat(:)
    complex(kind=dbl),       allocatable :: velc_cp(:)
    
    allocate( velc_cp(this%jmv), heat(this%nd) )
    
      do ir = 1, this%nd
        velc_cp = this%rad_grid%c(ir,-1) * this%cp_fn(ir  ) * this%sol%velocity_jml_fn(ir  ) + &
                & this%rad_grid%c(ir,+1) * this%cp_fn(ir+1) * this%sol%velocity_jml_fn(ir+1)
             
        heat(ir) = dotproduct_fn( this%jmax, velc_cp, this%sol%flux_jml_fn(ir) / this%lambda_fn(ir) )
      end do
      
    deallocate( velc_cp )
    
      advected_heat_fn = this%rad_grid%intV_fn( heat )
    
    deallocate( heat )
    
  end function advected_heat_fn
   
end submodule BalanceEquations