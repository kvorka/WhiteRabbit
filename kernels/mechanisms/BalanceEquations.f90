module BalanceEquations
  use PhysicalObject
  use NonLinearTerms
  implicit none
  
  public :: laws_mech_fn
  public :: laws_force_fn
  public :: laws_temp_fn
  
  private :: buoyancy_power_fn
  private :: heating_power_fn
  private :: bound_power_fn
  private :: bound_flux_fn
  private :: advected_heat_fn
  
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
        laws_temp_fn = real(this%sol%flux_fn(this%nd,0,0,1), kind=dbl) / real(this%sol%flux_fn(1,0,0,1), kind=dbl) / (this%r_ud**2)
    end select
    
  end function laws_temp_fn
  
  real(kind=dbl) function laws_force_fn(this, j, m)
    class(T_physicalObject),          intent(in) :: this
    integer,                          intent(in) :: j, m
    integer                                      :: i
    complex(kind=dbl)                            :: press_topo_d, press_topo_u, press_buoy
    complex(kind=dbl), dimension(:), allocatable :: drho
    
    
    associate( rd => this%rad_grid%r(1),                          &
             & ru => this%rad_grid%r(this%nd),                    &
             & gd => this%gravity%g_fn(this%rad_grid%r(1)),       &
             & gu => this%gravity%g_fn(this%rad_grid%r(this%nd))  )
    
    press_topo_u = this%Rau * ru**2 * gu * this%sol%t_up(jm(j,m))
    press_topo_d = this%Rad * rd**2 * gd * this%sol%t_dn(jm(j,m))
    
    allocate( drho(this%nd+1) )
      do i = 1, this%nd+1
        drho(i) = this%Ra * this%alpha_fn(i) * this%gravity%g_fn(this%rad_grid%rr(i)) * this%sol%temp_fn(i,j,m)
      end do
      
      press_buoy = this%rad_grid%intV_fn( drho )
    deallocate( drho )
    
    end associate
    
    laws_force_fn = real( press_topo_u / ( press_topo_d + press_buoy ), kind=dbl )
    
  end function laws_force_fn
  
  
  real(kind=dbl) function buoyancy_power_fn(this)
    class(T_physicalObject),       intent(in) :: this
    integer                                   :: i
    real(kind=dbl), dimension(:), allocatable :: power_i
    
    allocate( power_i(this%nd+1) )
    
      do i = 1, this%nd+1
        power_i(i) = dotproduct_fn( this%jmax, erT_fn(this,i), this%sol%velocity_jml_fn(i) )
      end do
      
      buoyancy_power_fn = this%Ra * this%rad_grid%intV_fn( power_i )
    
    deallocate( power_i )

  end function buoyancy_power_fn
  
  real(kind=dbl) function heating_power_fn(this)
    class(T_physicalObject),       intent(in) :: this
    integer                                   :: i
    real(kind=dbl), dimension(:), allocatable :: power_i
    
    allocate( power_i(this%nd) )
    
      do i = 1, this%nd
        power_i(i) = tnorm_fn(this%jmax, this%sol%deviatoric_stress_jml2_fn(i))**2 / this%visc_fn(i) / 2
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
    
    bound_flux_fn = ( real(-this%sol%flux_fn(this%nd,0,0,+1), kind=dbl) * this%rad_grid%r(this%nd)**2 - &
                    & real(-this%sol%flux_fn(1,0,0,+1), kind=dbl) * this%rad_grid%r(1)**2 - &
                    & this%Ds/this%Ra * this%rad_grid%intV_fn( real(this%htide(:,1),kind=dbl) ) ) / sqrt(4*pi)
    
  end function bound_flux_fn
  
  real(kind=dbl) function advected_heat_fn(this)
    class(T_physicalObject),          intent(in) :: this
    integer                                      :: i
    real(kind=dbl),    dimension(:), allocatable :: heat
    complex(kind=dbl), dimension(:), allocatable :: velc_cp
    
    allocate( velc_cp(this%jmv), heat(this%nd) )
    
      do i = 1, this%nd
        velc_cp = this%rad_grid%c(i,-1) * this%cp_fn(i  ) * this%sol%velocity_jml_fn(i  ) + &
                & this%rad_grid%c(i,+1) * this%cp_fn(i+1) * this%sol%velocity_jml_fn(i+1)
             
        heat(i) = dotproduct_fn( this%jmax, velc_cp, this%sol%flux_jml_fn(i) / this%lambda_fn(i) )
      end do
      
    deallocate( velc_cp )
    
      advected_heat_fn = this%rad_grid%intV_fn( heat )
    
    deallocate( heat )
    
  end function advected_heat_fn
   
end module BalanceEquations