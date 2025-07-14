module ocean_constants
  use math
  implicit none
  
  !Vseobecne nastavenie konstant
  integer, parameter :: nd_ocean = 81
  integer, parameter :: jmax_ocean = 237
  integer, parameter :: n_iter_ocean = 200  !Nastavit na N_period pre slapy
  logical, parameter :: noharm_ocean = .false.
  
  character(len=*), parameter :: grid_type_ocean    = 'chebv'
  character(len=*), parameter :: rheology_ocean     = 'viscos'
  character(len=*), parameter :: thermal_bnd_ocean  = 'fluxd'
  character(len=*), parameter :: mechanic_bnd_ocean = 'frees'
  character(len=*), parameter :: gravity_ocean      = 'new'
  character(len=*), parameter :: scaling_ocean      = 'basics'
  
  !Nastavovanie konstant pre konvektivny vypocet
  real(kind=dbl), parameter :: r_ud_ocean = 0.90_dbl
  real(kind=dbl), parameter :: D_ud_ocean = 1.0d5
  
  !Nastavovanie konstant pre slapovy vypocet
  real(kind=dbl), parameter, private :: rd         = 194100  !m
  real(kind=dbl), parameter, private :: ru         = 241100  !m
  real(kind=dbl), parameter          :: period     = 118387  !s
  real(kind=dbl), parameter, private :: nu         = 1e3     !Pa.s
  real(kind=dbl), parameter, private :: rho        = 1e3     !kg/m3
  real(kind=dbl), parameter          :: stress_dim = D_ud_ocean**3 * nu * (2*pi / period)**2
  
  real(kind=dbl), parameter :: Pr_ocean = 10._dbl
  real(kind=dbl), parameter :: Ra_ocean = 2.50d7
  real(kind=dbl), parameter :: Ek_ocean = 6.0d-4
  real(kind=dbl), parameter :: Kl_ocean = zero
  
  !Nastavovanie pociatocneho stavu
  logical, parameter :: init_through_file_ocean = .true.
  logical, parameter :: init_through_file_bnd_ocean = .false.
  integer, parameter :: nd_init_ocean = 81
  integer, parameter :: jmax_init_ocean = 237
  
end module ocean_constants
