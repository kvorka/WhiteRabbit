program ConvergenceIce
  use ConvergenceCurveMod
  implicit none
  
  call convergence_curve_ice_sub('curveShape.dat', 'code/ice/crust/data_shape/Shape-')
  call convergence_curve_ice_sub('curveTopo.dat',  'code/ice/crust/data_topo/Topo-')
  
end program ConvergenceIce
  