submodule (ice) init
  implicit none ; contains
  
  module procedure init_ice_sub
    
    if (present(noharm)) then
      call this%init_objects_sub( nd = nd_ice, jmax = jmax_in, r_ud = rdown_ice / rup_ice, rgrid = grid_type_ice, &
                                & gmod = gravity_ice, g = g_ice, noharm = noharm                                  )
    else
      call this%init_objects_sub( nd = nd_ice, jmax = jmax_in, r_ud = rdown_ice / rup_ice, rgrid = grid_type_ice, &
                                & gmod = gravity_ice, g = g_ice                                                   )
    end if
    
    this%n_iter = n_iter
    
    this%rheology     = rheol_in
    this%mechanic_bnd = mechanic_bnd_ice
    this%thermal_bnd  = thermal_bnd_ice

    this%Pr = huge(zero)
    this%Ek = huge(zero)

    this%g    = g_ice
    this%D_ud = D_ice
    this%Td   = name_meltingTemp_fn( D_ice )
    this%Tu   = theta_average_fn( name_surfaceTemp_fn, 0.0054_dbl )
    
    this%rhoI  = rho_ice
    this%rhoW  = rho_water
    this%rhoI2 = rho_iceII
    this%rhoC  = rho_core
    this%rC    = r_core  / this%D_ud
    this%rI2   = r_iceII / this%D_ud
    
    this%diam    = diam_ice
    this%cutoff  = cutoff_ice
    this%mu      = mu_ice
    this%lambdaC = lambdaC_ice
    this%hC      = hC_ice / this%D_ud
    this%omega   = omega
    
    this%alphaU  = name_expansivity_fn( this%Tu )
    this%lambdaU = name_conductivity_fn( this%Tu )
    this%cU      = name_cp_fn( this%Tu )
    this%viscU   = goldsby_visc_fn( this%diam, (this%Tu+this%Td)/2, zero )
    this%kappaU  = this%lambdaU / ( this%cU * this%rhoI )
    
    this%period = 2 * pi / omega * ( this%kappaU / this%D_ud**2 )
    
    this%Ra   = (this%rhoI * this%alphaU * (this%Td-this%Tu)) * this%g * this%D_ud**3 / ( this%viscU * this%kappaU )
    this%Raf  =                  this%cU * (this%Td-this%Tu) / lI_ice
    this%Ramu = this%viscU * this%kappaU / this%D_ud**2 / this%mu
    this%Rad  = (this%rhoI-this%rhoW) * this%g * this%D_ud**3 / ( this%viscU * this%kappaU )
    this%Rau  = (this%rhoI          ) * this%g * this%D_ud**3 / ( this%viscU * this%kappaU )
    this%Ds   = this%alphaU * this%g * this%D_ud / this%cU
    this%Cl   = this%g * this%D_ud * (this%rhoW-this%rhoI) * this%Td / ( this%rhoI * lI_ice * (this%Td-this%Tu) )
    
    this%gam = gamma_ice * this%D_ud**2 / this%lambdaU / (this%Td-this%Tu)
    
    call this%gravity%set_sub( Dcrust = this%D_ud, omega = omega, exc = exc )
    call this%bnd%init_layers_sub()
    
  end procedure init_ice_sub
  
  module procedure deallocate_ice_sub
    
    if ( allocated(this%nsph1) ) deallocate( this%nsph1 )
    if ( allocated(this%nsph2) ) deallocate( this%nsph2 )
    if ( allocated(this%ntorr) ) deallocate( this%ntorr )
    if ( allocated(this%ntemp) ) deallocate( this%ntemp )
    if ( allocated(this%nflux) ) deallocate( this%nflux )
    
    call this%gravity%deallocate_sub()
    call this%deallocate_objects_sub()
    
  end procedure deallocate_ice_sub
  
end submodule init
