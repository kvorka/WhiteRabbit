program ConvergenceIce
  use Paths
  use OutputIceMod
  implicit none
  
  call convergence_curve_ice_sub(path_ice_shape_dn, 'curveShape_dn.dat')
  call convergence_curve_ice_sub(path_ice_shape_up, 'curveShape_up.dat')
  call convergence_curve_ice_sub(path_ice_topo_dn , 'curveTopo_dn.dat' )
  call convergence_curve_ice_sub(path_ice_topo_up , 'curveTopo_up.dat' )
  
end program ConvergenceIce