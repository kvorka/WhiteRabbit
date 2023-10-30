module OutputIceMod
  use IceConstants
  use OutputMod
  implicit none
  
  integer, parameter, public :: jms = jmax_ice * (jmax_ice+1) / 2 + jmax_ice + 1
  
  public :: convergence_curve_ice_sub, harm_analysis_ice_sub, save_spectra_ice_sub
  
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
    integer                        :: ij, ijm, error
    complex(kind=dbl), allocatable :: spectra(:)
    real(kind=dbl),    allocatable :: spectra_j0(:), griddata_1d(:), griddata_2d(:,:)
    
    allocate( spectra(jms) )
    
    open(unit=1, file=path_in, status='old', action='read')
      do
        read(1,*,iostat=error) ijm, spectra(ijm) ; if (error /= 0) exit
      end do
    close(1)
    
    if (zon) then
      allocate( spectra_j0(0:jmax_ice) , griddata_1d(nth) )
        
        do ij = 0, jmax_ice
          spectra_j0(ij) = c2r_fn( spectra( jm(ij,0) ) )
        end do

        call harmsy_Pj0_sub(jmax_ice, spectra_j0, griddata_1d)
        call out_data_1d_sub(path_out, griddata_1d)
      
      deallocate( spectra_j0 , griddata_1d, spectra )
    else
      allocate( griddata_2d(2*nth,nth) )
        
        call harmsy_sub(jmax_ice, spectra, griddata_2d)
        call out_data_2d_sub(path_out, griddata_2d)
        
      deallocate( griddata_2d, spectra )
    end if
    
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
