submodule (PhysicalObject) BalanceEquations
  implicit none ; contains
  
  module real(kind=dbl) function laws_mech_fn(this)
    class(T_physicalObject), intent(in) :: this
    
    select case( this%mechanic_bnd )
      case( 'shape' )
        laws_mech_fn = ( this%bottombnd_power_fn()  + this%upperbnd_power_fn() ) / &
                     & ( this%viscdissip_power_fn() - this%buoyancy_power_fn() )
        
      case default  
        laws_mech_fn = this%buoyancy_power_fn() / this%viscdissip_power_fn()
        
    end select
    
  end function laws_mech_fn
  
  module real(kind=dbl) function laws_temp_fn(this)
    class(T_physicalObject), intent(in) :: this
    integer                             :: ir
    real(kind=dbl)                      :: fac1, fac2, fac3, flux_dn, flux_up, totheat, totheattide
    real(kind=dbl),         allocatable :: heat(:), heattide(:)
    
    flux_dn = c2r_fn( +this%sol%flux_fn(1,1,1)       * this%rd**2 )
    flux_up = c2r_fn( -this%sol%flux_fn(this%nd,1,1) * this%ru**2 )
    
    select case( this%thermal_bnd )
      case( 'phase' )
        allocate( heat(this%nd), heattide(this%nd) )
          
          do ir = 1, this%nd
            fac1 = this%rad_grid%c(ir,-1) * this%cp_fn(ir  )
            fac2 = this%rad_grid%c(ir,+1) * this%cp_fn(ir+1)
            fac3 = 1 / this%lambda_fn(ir)
            
            heat(ir) = dotproduct_fn( this%jmax , fac1 * this%sol%velocity_jml_fn(ir)  + &
                                                & fac2 * this%sol%velocity_jml_fn(ir+1), &
                                                & fac3 * this%sol%flux_jml_fn(ir)        )
            
            heattide(ir) = this%htide_fn(ir,1)
          end do
          
          totheat     = this%rad_grid%intV_fn( heat )
          totheattide = this%rad_grid%intV_fn( heattide )
          
        deallocate( heat, heattide )
        
        laws_temp_fn = totheat / ( flux_dn + flux_up + totheattide )
        
      case default
        
        laws_temp_fn = -flux_up / flux_dn
        
    end select
    
  end function laws_temp_fn
  
  module real(kind=dbl) function laws_force_fn(this, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ijm
    integer                             :: ir
    complex(kind=dbl)                   :: press_topo_d, press_topo_u, press_buoy
    complex(kind=dbl),      allocatable :: dbuoy(:)
    
    press_topo_d = this%Rad * this%rd**2 * this%gd * this%sol%t_dn(ijm)
    press_topo_u = this%Rau * this%ru**2 * this%gu * this%sol%t_up(ijm)
    
    allocate( dbuoy(this%nd+1) )
    
      do ir = 1, this%nd+1
        dbuoy(ir) = this%Ra * this%alpha_fn(ir) * this%gravity%g_fn(this%rad_grid%rr(ir)) * this%sol%temp_fn(ir,ijm)
      end do
      
      press_buoy = this%rad_grid%intV_fn( dbuoy )
    
    deallocate( dbuoy )
    
    laws_force_fn = c2r_fn( press_topo_u / ( press_topo_d + press_buoy ) )
    
  end function laws_force_fn
  
end submodule BalanceEquations