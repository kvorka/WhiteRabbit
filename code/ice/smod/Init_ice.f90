submodule(IceMod) Init_ice
  implicit none ; contains
  
  module subroutine init_ice_sub(this, jmax_in, rheol_in, n_iter, noharm)
    class(T_ice),      intent(inout) :: this
    integer,           intent(in)    :: jmax_in, n_iter
    character(len=*),  intent(in)    :: rheol_in
    logical, optional, intent(in)    :: noharm
    
    if (present(noharm)) then
      call this%init_objects_sub( nd = nd_ice, jmax = jmax_in, r_ud = rdown_ice / rup_ice, rgrid = grid_type_ice, &
                                & gmod = gravity_ice, g = g_ice, noharm = noharm                                  )
    else
      call this%init_objects_sub( nd = nd_ice, jmax = jmax_in, r_ud = rdown_ice / rup_ice, rgrid = grid_type_ice, &
                                & gmod = gravity_ice, g = g_ice                                                   )
    end if
    
    this%n_iter = n_iter
    
    this%rheology     = rheol_in
    this%andrade      = .false.
    this%mechanic_bnd = mechanic_bnd_ice
    this%thermal_bnd  = thermal_bnd_ice

    this%Pr = huge(0._dbl)
    this%Ek = huge(0._dbl)

    this%g    = g_ice
    this%Td   = Td_ice
    this%Tu   = Tu_ice
    this%D_ud = rup_ice - rdown_ice
    
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
    
    this%alphaU  = 1.0d-4
    this%lambdaU = 0.4685_dbl + 488.12_dbl / this%Tu
    this%cU      = 185._dbl + 7.037_dbl * this%Tu
    this%viscU   = (this%Tu+this%Td)/2 * this%diam**2 * exp( 59.0d+3 / rgas / ( (this%Tu+this%Td)/2 ) ) / 9.0d-8 / 2
    this%kappaU  = this%lambdaU / ( this%cU * this%rhoI )
    
    this%period = 2 * pi / omega * ( this%kappaU / this%D_ud**2 )
    
    this%Ra   = (this%rhoI * this%alphaU * (this%Td-this%Tu)) * this%g * this%D_ud**3 / ( this%viscU * this%kappaU )
    this%Raf  =                  this%cU * (this%Td-this%Tu) / lI_ice
    this%Ramu = this%viscU * this%kappaU / this%D_ud**2 / this%mu
    this%Rad  = (this%rhoI-this%rhoW) * this%g * this%D_ud**3 / ( this%viscU * this%kappaU )
    this%Rau  = (this%rhoI          ) * this%g * this%D_ud**3 / ( this%viscU * this%kappaU )
    this%Ds   = this%alphaU * this%g * this%D_ud / this%cU
    this%Cl   = this%g * this%D_ud * (this%rhoW-this%rhoI) * this%Td / ( this%rhoI * lI_ice * (this%Td-this%Tu) )
    
    
    call this%gravity%set_sub( Dcrust = this%D_ud, omega = omega, exc = exc )
    call this%sol%init_layers_sub()
    
    allocate( this%htide(this%nd,jms4) ); this%htide = czero

  end subroutine init_ice_sub
  
  module subroutine deallocate_ice_sub(this)
    class(T_ice), intent(inout) :: this
    
    call this%gravity%deallocate_sub()
    call this%deallocate_objects_sub()
    
  end subroutine deallocate_ice_sub
  
end submodule Init_ice
