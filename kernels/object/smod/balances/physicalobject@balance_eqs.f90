submodule (physicalobject) balance_eqs
  implicit none ; contains
  
  module procedure laws_mech_fn
    integer                        :: ir, ijm
    real(kind=dbl)                 :: bndpow, heatpow, buoypow
    real(kind=dbl),    allocatable :: power_i(:)
    complex(kind=dbl), allocatable :: gdrho_jm(:), rvelc_jm(:), devstress_jm(:)
    
    !Viscous dissipation
    allocate( power_i(this%nd), devstress_jm(this%jmt) )
    
      !$omp parallel do private (devstress_jm)
      do ir = 1, this%nd
        devstress_jm = this%sol%deviatoric_stress_jml2_fn(ir)
        power_i(ir)  = tensnorm2_fn( this%jmax, devstress_jm ) / this%visc_r_fn(ir) / 2
      end do
      !$omp end parallel do
      
      heatpow = this%rad_grid%intV_fn( power_i )
    
    deallocate( power_i )
    
    !Buoyancy power
    allocate( power_i(this%nd+1) , gdrho_jm(this%jms), rvelc_jm(this%jms) )
    
      !$omp parallel do private (gdrho_jm, rvelc_jm)
      do ir = 1, this%nd+1
        call this%er_buoy_rr_jm_sub(ir, gdrho_jm)
        call this%vr_rr_jm_sub(ir, rvelc_jm)
        
        power_i(ir) = scalproduct_fn( this%jmax, gdrho_jm, rvelc_jm )
      end do
      !$omp end parallel do
      
      buoypow = this%rad_grid%intV_fn( power_i )
    
    deallocate( power_i, gdrho_jm )
    
    select case( this%mechanic_bnd )
      case( 'shape' )
        !Power of the bottom boundary
        do concurrent ( ijm = 1:this%jms )
          rvelc_jm(ijm) = this%vr_r_fn(1,ijm)
        end do
        
        bndpow = this%Rad * this%gd * this%rd**2 * scalproduct_fn(this%jmax, this%bnd%t_dn, rvelc_jm)
        
        !Power of the upper boundary
        do concurrent ( ijm = 1:this%jms )
          rvelc_jm(ijm) = this%vr_r_fn(this%nd,ijm)
        end do
        
        bndpow = bndpow - this%Rau * this%gu * this%ru**2 * scalproduct_fn(this%jmax, this%bnd%t_up, rvelc_jm)
        
        !Resulting law
        laws_mech_fn = bndpow / ( heatpow - buoypow )
        
      case default  
        laws_mech_fn = buoypow / heatpow
        
    end select
    
    deallocate( rvelc_jm )
    
  end procedure laws_mech_fn
  
  !!TODO: implement non-constant capacity case
  module procedure laws_temp_fn
    integer                        :: ir
    real(kind=dbl)                 :: flow_dn, flow_up, totheat, totheattide
    real(kind=dbl),    allocatable :: heat(:), heattide(:)
    complex(kind=dbl), allocatable :: velocity(:), gradT(:)
    
    flow_dn = c2r_fn( +this%q_r_fn(      1,1,1) ) * this%rd**2
    flow_up = c2r_fn( -this%q_r_fn(this%nd,1,1) ) * this%ru**2
    
    select case( this%thermal_bnd )
      case( 'phase' )
        allocate( heat(this%nd), heattide(this%nd), velocity(this%jmv), gradT(this%jmv) )
          
          !$omp parallel do private (velocity, gradT)
          do ir = 1, this%nd
            call this%v_r_ijml_sub( ir, velocity )
            call this%gradT_r_ijml_sub( ir, gradT, -1 )
            
            heat(ir)     = dotproduct_fn( this%jmax , velocity , gradT )
            heattide(ir) = c2r_fn( this%htide_r_fn(ir,1) )
          end do
          !$omp end parallel do
          
          totheat     = this%rad_grid%intV_fn( heat )
          totheattide = this%rad_grid%intV_fn( heattide )
          
        deallocate( heat, heattide, velocity, gradT )
        
        laws_temp_fn = totheat / ( flow_dn + flow_up + totheattide )
        
      case default
        
        laws_temp_fn = -flow_up / flow_dn
        
    end select
    
  end procedure laws_temp_fn
  
  module procedure laws_force_fn
    integer                        :: ir
    complex(kind=dbl)              :: press_topo_d, press_topo_u, press_buoy
    complex(kind=dbl), allocatable :: dbuoy(:)
    
    press_topo_d = this%Rad * this%rd**2 * this%gd * this%bnd%t_dn(ijm)
    press_topo_u = this%Rau * this%ru**2 * this%gu * this%bnd%t_up(ijm)
    
    allocate( dbuoy(this%nd+1) )
    
      do ir = 1, this%nd+1
        dbuoy(ir) = this%Ra * this%alpha_rr_fn(ir) * this%gravity%g_fn(this%rad_grid%rr(ir)) * this%temp_rr_fn(ir,ijm)
      end do
      
      press_buoy = this%rad_grid%intV_fn( dbuoy )
    
    deallocate( dbuoy )
    
    laws_force_fn = c2r_fn( press_topo_u / ( press_topo_d + press_buoy ) )
    
  end procedure laws_force_fn
  
end submodule balance_eqs