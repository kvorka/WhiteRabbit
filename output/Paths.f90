module Paths
  use Math
  implicit none
  
  real(kind=dbl), parameter :: tNuss      = 0.01_dbl
  integer,        parameter :: avrg_start = 5600
  integer,        parameter :: avrg_end   = 6600
  
  character(len=*), parameter :: path_nuss         = 'data/Nuss.dat'
  character(len=*), parameter :: path_ocean_temp   = 'data/data_ocean_temp/Temp-'
  character(len=*), parameter :: path_ocean_velc   = 'data/data_ocean_veloc/Velc-'
  character(len=*), parameter :: path_ocean_flux   = 'data/data_ocean_flux/Flux-'
  character(len=*), parameter :: path_ice_shape_dn = 'data/data_ice_shape/Shape_dn-'
  character(len=*), parameter :: path_ice_topo_dn  = 'data/data_ice_topo/Topo_dn-'
  character(len=*), parameter :: path_ice_shape_up = 'data/data_ice_shape/Shape_up-'
  character(len=*), parameter :: path_ice_topo_up  = 'data/data_ice_topo/Topo_up-'

end module Paths