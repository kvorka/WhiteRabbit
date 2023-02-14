module iceMod
  use Math
  use PhysicalObject
  use IceConstants
  implicit none
  
  type, extends(T_physicalObject), abstract, public :: T_ice
    real(kind=dbl) :: diam, lambdaC, hC, lambdaU, viscU, cutoff, alphaU, cU, Td, Tu, period, g
    real(kind=dbl) :: rC, rI2, rhoC, rhoI2, rhoW, rhoI
    
    contains
    
    procedure :: init_ice_sub       => init_ice_sub
    procedure :: lambda_fn          => lambda_ice_fn
    procedure :: cp_fn              => cp_ice_fn
    procedure :: alpha_fn           => alpha_ice_fn
    procedure :: visc_fn            => visc_ice_fn
    procedure :: deallocate_ice_sub => deallocate_ice_sub
    
  end type T_ice
  
  private :: init_ice_sub
  private :: deallocate_ice_sub
  
  private :: lambda_ice_fn
  private :: cp_ice_fn
  private :: alpha_ice_fn
  private :: visc_ice_fn
  
  contains

  subroutine init_ice_sub(this, jmax_in, rheol_in, n_iter, noharm)
    class(T_ice),      intent(inout) :: this
    integer,           intent(in)    :: jmax_in, n_iter
    character(len=*),  intent(in)    :: rheol_in
    logical, optional, intent(in)    :: noharm
    real(kind=dbl)                   :: tkappa_ice
    
    !Inicializuj objekty vseobecne
    if (present(noharm)) then
      call this%init_objects_sub( nd = nd_ice, jmax = jmax_in, r_ud = rdown_ice / rup_ice, rgrid = grid_type_ice, noharm = noharm )
    else
      call this%init_objects_sub( nd = nd_ice, jmax = jmax_in, r_ud = rdown_ice / rup_ice, rgrid = grid_type_ice )
    end if
    
    this%n_iter = n_iter
    
    this%rheology     = rheol_in
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
    this%lambdaC = lambdaC_ice
    this%hC      = hC_ice / this%D_ud
    
    this%alphaU  = 1.0d-4
    this%lambdaU = 0.4685_dbl + 488.12_dbl / this%Tu
    this%cU      = 185._dbl + 7.037_dbl * this%Tu
    this%viscU   = (this%Tu+this%Td)/2 * this%diam**2 * exp( 59.0d+3 / rgas / ( (this%Tu+this%Td)/2 ) ) / 9.0d-8 / 2
    
    tkappa_ice = this%lambdaU / this%cU / this%rhoI
    
    this%period = 2 * pi / omega * ( tkappa_ice / this%D_ud**2 )
    
    this%Ra   = (this%rhoI * this%alphaU * (this%Td-this%Tu)) * this%g * this%D_ud**3 / this%viscU / tkappa_ice
    this%Raf  =                  this%cU * (this%Td-this%Tu) / lI_ice
    this%Ramu = mu_ice * this%D_ud**2                                                 / this%viscU / tkappa_ice
    this%Rad  = (this%rhoI-this%rhoW                        ) * this%g * this%D_ud**3 / this%viscU / tkappa_ice
    this%Rau  = (this%rhoI                                  ) * this%g * this%D_ud**3 / this%viscU / tkappa_ice
    this%Ds   = this%alphaU * this%g * this%D_ud / this%cU
    
    !Inicializuj gravitaciu
    call this%gravity%init_sub( gmod = gravity_ice, g = this%g, Dcrust = this%D_ud, omega = omega, exc = exc )
    
    !Inicializuj premenne pre riesenie
    call this%sol%init_stemp_sub()
    call this%sol%init_smech_sub()
    call this%sol%init_layers_sub()
    
    !Inicializuj premenne pre matice
    call this%mat%init_mtemp_sub()
    call this%mat%init_mmech_sub()
    
    !Inicializuj slapove zahrievanie
    allocate( this%htide(this%nd,jms4) ); this%htide = cmplx(0._dbl, 0._dbl, kind=dbl)
    
  end subroutine init_ice_sub
  
  real(kind=dbl) function lambda_ice_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    real(kind=dbl)            :: temp, lambdaI

    if ( this%rad_grid%r(i) < this%rad_grid%r(this%nd) - this%hC ) then
       temp = this%Tu + (this%Td-this%Tu) * real( this%rad_grid%c(i,-1) * this%sol%temp_fn(i  ,0,0) + &
                                                & this%rad_grid%c(i,+1) * this%sol%temp_fn(i+1,0,0)    , kind=dbl) / sqrt(4*pi)
       
       lambdaI = 0.4685_dbl + 488.12_dbl / temp
       
    else
       lambdaI = this%lambdaC
       
    end if
     
    lambda_ice_fn = lambdaI / this%lambdaU
    
  end function lambda_ice_fn
  
  pure real(kind=dbl) function cp_ice_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    real(kind=dbl)            :: temp, cI
    
    temp = this%Tu + (this%Td-this%Tu) * real( this%sol%temp_fn(i,0,0), kind=dbl ) / sqrt(4*pi)
    cI = 185._dbl + 7.037_dbl * temp
    
    cp_ice_fn = cI / this%cU
    
  end function cp_ice_fn
  
  pure real(kind=dbl) function alpha_ice_fn(this, i)
    real(kind=dbl), parameter :: A0 = 128.2147_dbl
    real(kind=dbl), parameter :: A1 = 0._dbl
    real(kind=dbl), parameter :: A2 = 0._dbl
    real(kind=dbl), parameter :: A3 = -1.3152d-6
    real(kind=dbl), parameter :: A4 = +2.4837d-8
    real(kind=dbl), parameter :: A5 = -1.6064d-10
    real(kind=dbl), parameter :: A6 = +4.6097d-13
    real(kind=dbl), parameter :: A7 = -4.9661d-16
    
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    real(kind=dbl)            :: temp, v, dv, alphaI
    
    temp = this%Tu + (this%Td-this%Tu) * real( this%sol%temp_fn(i,0,0), kind=dbl ) / sqrt(4*pi)
    
    v  = A0 +    A3 * temp**3 +     A4 * temp**4 +     A5 * temp**5 +     A6 * temp**6 +     A7 * temp**7
    dv =      3* A3 * temp**2 + 4 * A4 * temp**3 + 5 * A5 * temp**4 + 6 * A6 * temp**5 + 7 * A7 * temp**6
    
    alphaI = dv/v
    
    alpha_ice_fn = alphaI / this%alphaU
    
  end function alpha_ice_fn
  
  pure function visc_ice_fn(this, i) result(viscI)
    real(kind=dbl), parameter :: a = 9.0d-8
    real(kind=dbl), parameter :: e = 59.0d+3
    
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    real(kind=dbl)            :: temp, viscI
    
    temp = this%Tu + (this%Td-this%Tu) * real( this%rad_grid%c(i,-1) * this%sol%temp_fn(i,0,0) + &
                                             & this%rad_grid%c(i,+1) * this%sol%temp_fn(i+1,0,0), kind=dbl) / sqrt(4*pi)
    
    viscI = min( temp * this%diam**2 * exp(e / rgas / temp) / a / 2, this%cutoff )
    
    viscI = viscI / this%viscU
    
  end function visc_ice_fn

  subroutine deallocate_ice_sub(this)
    class(T_ice), intent(inout) :: this
    
    deallocate( this%htide )
    
    call this%deallocate_objects_sub()
    
  end subroutine deallocate_ice_sub
  
end module iceMod
