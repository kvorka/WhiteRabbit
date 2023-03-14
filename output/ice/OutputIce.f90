program OutputIce
  use Paths
  use OutputIceMod
  implicit none

  call save_spectra_ice_sub(path_ice_shape, 'Shape')
  call save_spectra_ice_sub(path_ice_topo , 'Topo' )
  
  call harm_analysis_ice_sub('Shape')
  call harm_analysis_ice_sub('Topo')

  call zonal_analysis_ice_sub('Shape')
  call zonal_analysis_ice_sub('Topo')

end program OutputIce
