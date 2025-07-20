module Paths
  use math
  implicit none
  
  real(kind=dbl), parameter :: tNuss      = 0.26_dbl
  integer,        parameter :: avrg_start = 1300
  integer,        parameter :: avrg_end   = 1750
  
  character(len=*), parameter :: path_nuss         = 'data/Nuss.dat'
  character(len=*), parameter :: path_ocean_temp   = 'data/data_ocean_temp/Temp-'
  character(len=*), parameter :: path_ocean_velc   = 'data/data_ocean_veloc/Velc-'
  character(len=*), parameter :: path_ocean_flux   = 'data/data_ocean_fluxu/Fluxu-'
  character(len=*), parameter :: path_ice_shape_dn = 'data/data_ice_shape/Shape_dn-'
  character(len=*), parameter :: path_ice_topo_dn  = 'data/data_ice_topo/Topo_dn-'
  character(len=*), parameter :: path_ice_shape_up = 'data/data_ice_shape/Shape_up-'
  character(len=*), parameter :: path_ice_topo_up  = 'data/data_ice_topo/Topo_up-'

end module Paths
