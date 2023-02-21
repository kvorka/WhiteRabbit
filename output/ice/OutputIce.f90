program OutputIce
  use OutputIceMod
  implicit none

  call save_spectra_ice_sub('code/ice/crust/data_shape/Shape-', 'Shape')
  call save_spectra_ice_sub('code/ice/crust/data_shape/Topo-', 'Topo')
  
  call harm_analysis_ice_sub('Shape')
  call harm_analysis_ice_sub('Topo')

end program OutputIce
