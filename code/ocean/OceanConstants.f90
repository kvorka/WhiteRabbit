module OceanConstants
    use Math
    implicit none
    
    !Vseobecne nastavenie konstant
    integer, parameter :: nd_ocean = 73
    integer, parameter :: jmax_ocean = 125
    integer, parameter :: n_iter_ocean = 200  !Nastavit na N_period pre slapy
    logical, parameter :: noharm_ocean = .false.

    character(len=*), parameter :: grid_type_ocean    = 'chebv'
    character(len=*), parameter :: rheology_ocean     = 'viscos'
    character(len=*), parameter :: thermal_bnd_ocean  = 'basic'
    character(len=*), parameter :: mechanic_bnd_ocean = 'noslp'
    character(len=*), parameter :: gravity_ocean      = 'new'
    character(len=*), parameter :: scaling_ocean      = 'basics'

    !Nastavovanie konstant pre konvektivny vypocet
    real(kind=dbl), parameter :: r_ud_ocean = 0.6_dbl
    real(kind=dbl), parameter :: D_ud_ocean = 1.0d5

    !Nastavovanie konstant pre slapovy vypocet
    real(kind=dbl), parameter, private :: rd         = 194100  !m
    real(kind=dbl), parameter, private :: ru         = 241100  !m
    real(kind=dbl), parameter          :: period     = 118387  !s
    real(kind=dbl), parameter, private :: nu         = 1e3     !Pa.s
    real(kind=dbl), parameter, private :: rho        = 1e3     !kg/m3
    real(kind=dbl), parameter          :: stress_dim = D_ud_ocean**3 * nu * (2*pi / period)**2
    
    real(kind=dbl), parameter :: Pr_ocean = 1._dbl
    real(kind=dbl), parameter :: Ra_ocean = 1.0d7
    real(kind=dbl), parameter :: Ek_ocean = 1.0d-4
    !real(kind=dbl), parameter :: Ek_ocean = nu / rho * period / 2 / pi / D_ud_ocean**2
    real(kind=dbl), parameter :: Cl_ocean = 0._dbl
    
    !Nastavovanie pociatocneho stavu
    logical, parameter :: init_through_file_ocean = .false.
    integer, parameter :: nd_init_ocean = 73
    integer, parameter :: jmax_init_ocean = 186
    
  end module OceanConstants
