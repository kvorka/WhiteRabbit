module OutputIceMod
  use IceConstants
  use OutputMod
  implicit none
  
  integer, parameter, public :: jms = jmax_ice * (jmax_ice+1) / 2 + jmax_ice + 1
  
  public :: convergence_curve_ice_sub
  public :: harm_analysis_ice_sub
  public :: save_spectra_ice_sub
  
  contains
  
  subroutine convergence_curve_ice_sub(path_in, path_out)
    character(len=*), intent(in) :: path_out, path_in
    integer                      :: in, ij, im, error
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
              read(7,*,iostat=error) ij, im, coeff ; if (error /= 0) exit
              
              if (im == 0) then
                L2 = L2 + abs( coeff )**2
              else
                L2 = L2 + 2 * abs( coeff )**2
              end if
            end do
          end if
          
          write(8,*) in, sqrt(L2)
        close(7)
        
        in = in + 1
      end do
    close(8)
    
  end subroutine convergence_curve_ice_sub
  
  subroutine harm_analysis_ice_sub(path_in, path_out, zon)
    character(len=*),   intent(in) :: path_in, path_out
    logical,            intent(in) :: zon
    integer                        :: ijm, error
    complex(kind=dbl), allocatable :: spectra(:)
    real(kind=dbl),    allocatable :: griddata(:,:)
    
    allocate( spectra(jms) , griddata(2*nth,nth) )
    
    open(unit=1, file=path_in, status='old', action='read')
      do
        read(1,*,iostat=error) ijm, spectra(ijm) ; if (error /= 0) exit
      end do
    close(1)
    
    call harmsy_sub(jmax_ice, spectra, griddata)
    
    if (zon) then
      call out_data_1d_sub(path_out, griddata)
    else
      call out_data_2d_sub(path_out, griddata)
    end if
    
    deallocate( spectra, griddata )
    
  end subroutine harm_analysis_ice_sub

  subroutine save_spectra_ice_sub(path_in, path_out)
    character(len=*),  intent(in)  :: path_in, path_out
    complex(kind=dbl), allocatable :: spectra(:)
    
    allocate( spectra(jms) ); spectra = czero

    call avrg_spectra_2d_sub(path_in, spectra)
    call out_spectra_2d_sub(path_out, spectra)
    
    deallocate( spectra )
    
  end subroutine save_spectra_ice_sub
    
end module OutputIceMod
