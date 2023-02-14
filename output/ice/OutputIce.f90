program OutputIce
  use OutputIceMod
  implicit none
  
  call harm_analysis_ice_sub('code/ice/crust/data_shape/Shape-', 'Shape')
  call harm_analysis_ice_sub('code/ice/crust/data_topo/Topo-', 'Topo')

end program OutputIce
