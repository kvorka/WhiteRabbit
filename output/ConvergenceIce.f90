program ConvergenceIce
  use IceConstants
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
    real(kind=dbl)               :: L2
    complex(kind=dbl)            :: coeff
    
    open(unit=8, file=path_out, status='new', action='write'); in = 0
      do
        open(unit=7, file=path_in//trim(adjustl(int2str_fn(in)))//'.dat', status='old', action='read', iostat=error)
        
          if (error /= 0) then 
            exit
          else
            L2 = 0._dbl
            
            do
              read(7,*,iostat=error) ijm, coeff ; if (error /= 0) exit
              
              L2 = L2 + abs( coeff )**2
            end do
          end if
          
          write(8,*) in, sqrt(L2)
        close(7)
        
        in = in + 1
      end do
    close(8)
    
  end subroutine convergence_curve_ice_sub
  
end program ConvergenceIce