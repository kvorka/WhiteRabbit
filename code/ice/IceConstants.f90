module IceConstants
  use Math
  implicit none
  
  character(len=5), parameter :: grid_type_ice    = 'homog'
  character(len=5), parameter :: mechanic_bnd_ice = 'shape'
  character(len=5), parameter :: thermal_bnd_ice  = 'phase'
  character(len=3), parameter :: gravity_ice      = 'hom'
  
  integer, parameter :: nd_ice = 97
  integer, parameter :: jmax_ice = 24
  integer, parameter :: n_iter_ice = 10
  integer, parameter :: n_iter_tides = 1000
  
  integer, parameter :: nlay = 13
  integer, parameter :: jms4 = 15
  
  real(kind=dbl), parameter :: r_core  = 2068160._dbl
  real(kind=dbl), parameter :: r_iceII = 2149160._dbl
  real(kind=dbl), parameter :: rdown_ice = 2486040._dbl
  real(kind=dbl), parameter :: rup_ice   = 2574700._dbl
  
  real(kind=dbl), parameter :: rho_ice   = 934._dbl
  real(kind=dbl), parameter :: rho_water = 1200._dbl
  real(kind=dbl), parameter :: rho_iceII = 1371._dbl  
  real(kind=dbl), parameter :: rho_core  = 2532._dbl
  
  real(kind=dbl), parameter :: Td_ice = 265._dbl
  real(kind=dbl), parameter :: Tu_ice = 90._dbl
  
  real(kind=dbl), parameter :: g_ice      = 1.40_dbl
  real(kind=dbl), parameter :: diam_ice   = 1.0d-2
  real(kind=dbl), parameter :: cutoff_ice = 1.0d24
  
  real(kind=dbl), parameter :: lI_ice = 330000._dbl
  real(kind=dbl), parameter :: hC_ice = 5.0d3
  real(kind=dbl), parameter :: lambdaC_ice = 0.5_dbl
  real(kind=dbl), parameter :: mu_ice      = 3.5d9
  
  real(kind=dbl), parameter :: omega  = 4.5608d-6
  real(kind=dbl), parameter :: exc    = 0.02888_dbl
  
end module IceConstants
