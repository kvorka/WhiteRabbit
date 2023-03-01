program ConvergenceIce
  use Paths
  use ConvergenceCurveMod
  implicit none
  
  call convergence_curve_ice_sub('curveShape.dat', path_ice_shape)
  call convergence_curve_ice_sub('curveTopo.dat' , path_ice_topo )
  
end program ConvergenceIce