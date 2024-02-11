module OutputIceMod
  use IceConstants
  use OutputMod
  implicit none; public
  
  integer, parameter, private :: jms = jmax_ice * (jmax_ice+1) / 2 + jmax_ice + 1
  
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
  
  subroutine save_spectra_ice_sub(path_in, path_out)
    character(len=*),  intent(in)  :: path_in, path_out
    complex(kind=dbl), allocatable :: spectra(:)
    
    allocate( spectra(jms) ); spectra = czero

    call avrg_spectra_2d_sub(path_in, spectra)
    call out_spectra_2d_sub(path_out, spectra)
    
    deallocate( spectra )
    
  end subroutine save_spectra_ice_sub
    
end module OutputIceMod
