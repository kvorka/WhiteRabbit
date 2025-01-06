module ice_constants
  use math
  implicit none
  
  character(len=5), parameter :: grid_type_ice    = 'homog'
  character(len=5), parameter :: mechanic_bnd_ice = 'shape'
  character(len=5), parameter :: thermal_bnd_ice  = 'phase'
  character(len=3), parameter :: gravity_ice      = 'hom'
  
  integer, parameter :: nd_ice = 97
  integer, parameter :: jmax_ice = 47
  integer, parameter :: n_iter_ice = 10000
  integer, parameter :: n_iter_tides = 1000
  
  integer, parameter :: nlay = 13
  integer, parameter :: jms4 = 15
  
  real(kind=dbl), parameter :: D_ice     = 70.0d3
  real(kind=dbl), parameter :: rup_ice   = 2634.0d3
  real(kind=dbl), parameter :: rdown_ice = rup_ice - D_ice
  
  real(kind=dbl), parameter :: r_core  = 2068160._dbl
  real(kind=dbl), parameter :: r_iceII = 2149160._dbl
  
  real(kind=dbl), parameter :: rho_ice   = 920._dbl
  real(kind=dbl), parameter :: rho_water = 1200._dbl
  real(kind=dbl), parameter :: rho_iceII = 1371._dbl  
  real(kind=dbl), parameter :: rho_core  = 2532._dbl
  
  real(kind=dbl), parameter :: g_ice      = 1.43_dbl
  real(kind=dbl), parameter :: diam_ice   = 1.0d-2
  real(kind=dbl), parameter :: cutoff_ice = 1.0d24
  
  real(kind=dbl), parameter :: lI_ice = 334000._dbl
  real(kind=dbl), parameter :: hC_ice = 0._dbl
  real(kind=dbl), parameter :: lambdaC_ice = 0.5_dbl
  real(kind=dbl), parameter :: mu_ice      = 3.5d9
  
  real(kind=dbl), parameter :: omega  = 1.0164d-5
  real(kind=dbl), parameter :: exc    = 0.0015_dbl
  real(kind=dbl), parameter :: obl    = 0.033_dbl
  
  real(kind=dbl), parameter :: sun_to_body_dist = 7.780d11
  real(kind=dbl), parameter :: albedo = 0.43_dbl
  real(kind=dbl), parameter :: gamma_ice = 0._dbl
  
end module ice_constants
