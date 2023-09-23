program OutputIce
  use Paths
  use OutputIceMod
  implicit none

  call save_spectra_ice_sub(path_ice_shape_dn, 'shape_dn.spec')
  call save_spectra_ice_sub(path_ice_shape_up, 'shape_up.spec')
  call save_spectra_ice_sub(path_ice_topo_dn , 'topo_dn.spec' )
  call save_spectra_ice_sub(path_ice_topo_up , 'topo_up.spec' )
  
  call harm_analysis_ice_sub('shape_dn.spec', 'shape_dn.dat', zon=.false.)
  call harm_analysis_ice_sub('shape_up.spec', 'shape_up.dat', zon=.false.)
  call harm_analysis_ice_sub('topo_dn.spec' , 'topo_dn.dat' , zon=.false.)
  call harm_analysis_ice_sub('topo_up.spec' , 'topo_up.dat' , zon=.false.)

  call harm_analysis_ice_sub('shape_dn.spec', 'zon_shape_dn.dat', zon=.true.)
  call harm_analysis_ice_sub('shape_up.spec', 'zon_shape_up.dat', zon=.true.)
  call harm_analysis_ice_sub('topo_dn.spec' , 'zon_topo_dn.dat' , zon=.true.)
  call harm_analysis_ice_sub('topo_up.spec' , 'zon_topo_up.dat' , zon=.true.)

end program OutputIce
