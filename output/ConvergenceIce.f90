program ConvergenceIce
  use IceConstants
  use sph_norms
  use OutputMod
  implicit none
  
  call convergence_curve_ice_sub(path_ice_shape_dn, 'curveShape_dn.dat')
  call convergence_curve_ice_sub(path_ice_shape_up, 'curveShape_up.dat')
  call convergence_curve_ice_sub(path_ice_topo_dn , 'curveTopo_dn.dat' )
  call convergence_curve_ice_sub(path_ice_topo_up , 'curveTopo_up.dat' )
  
  contains
  
  subroutine convergence_curve_ice_sub(path_in, path_out)
    character(len=*), intent(in) :: path_out, path_in
    integer                      :: in, ijm, error
    complex(kind=dbl)            :: coeffs(jmax_ice*(jmax_ice+1)/2+jmax_ice+1)
    
    open(unit=8, file=path_out, status='new', action='write'); in = 0
    
    do
      open(unit=7, file=path_in//trim(adjustl(int2str_fn(in)))//'.dat', status='old', action='read', iostat=error)
      
        if (error /= 0) exit
        
        do ijm = 1, jmax_ice*(jmax_ice+1)/2+jmax_ice+1
          read(7,*) error, coeffs(ijm)
        end do
        
        write(8,*) in,  scalproduct_fn( jmax_ice, coeffs, coeffs )
      close(7)
      
      in = in + 1
    end do
    
    close(8)
    
  end subroutine convergence_curve_ice_sub
  
end program ConvergenceIce