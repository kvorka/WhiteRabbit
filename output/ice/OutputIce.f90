program OutputIce
  use Paths
  use OutputIceMod
  implicit none
  
  call save_spectra_ice_sub(path_ice_shape_dn, 'shape_dn.spec')
  call save_spectra_ice_sub(path_ice_shape_up, 'shape_up.spec')
  call save_spectra_ice_sub(path_ice_topo_dn , 'topo_dn.spec' )
  call save_spectra_ice_sub(path_ice_topo_up , 'topo_up.spec' )
  
end program OutputIce
