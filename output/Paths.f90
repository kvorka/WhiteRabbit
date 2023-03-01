module Paths
  implicit none

  integer, parameter :: avrg_start = 5600
  integer, parameter :: avrg_end   = 6600

  character(len=*), parameter :: path_nuss = 'data/Nuss.dat'

  character(len=*), parameter :: path_ocean_temp = 'data/data_ocean_temp/Temp-'
  character(len=*), parameter :: path_ocean_velc = 'data/data_ocean_veloc/Velc-'
  character(len=*), parameter :: path_ocean_flux = 'data/data_ocean_flux/Flux-'

  character(len=*), parameter :: path_ice_shape = 'data/data_ice_shape/Shape-'
  character(len=*), parameter :: path_ice_topo  = 'data/data_ice_topo/Topo-'

end module Paths