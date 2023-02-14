module OutputIceMod
  use Math
  use IceConstants
  implicit none

  complex(kind=dbl), allocatable :: spectra_up(:), spectra_dn(:)
  real(kind=dbl),    allocatable :: data_up(:,:), data_dn(:,:)
  
  public  :: harm_analysis_ice_sub
  private :: get_spectra_sub
 
  contains

  subroutine harm_analysis_ice_sub(path, opt)
    character(len=*), intent(in) :: path, opt
    
    allocate( spectra_up(jmax_ice*(jmax_ice+1)/2+jmax_ice+1) ); spectra_up = cmplx(0._dbl, 0._dbl, kind=dbl)
    allocate( spectra_dn(jmax_ice*(jmax_ice+1)/2+jmax_ice+1) ); spectra_dn = cmplx(0._dbl, 0._dbl, kind=dbl)
    
    call get_spectra_sub(path)
    
    allocate( data_up(2*nth,nth) ); data_up = 0._dbl
    allocate( data_dn(2*nth,nth) ); data_dn = 0._dbl
    
      !Preved spektrum na grid
      call toGrid_sub(jmax_ice, spectra_dn, data_dn)
      call toGrid_sub(jmax_ice, spectra_up, data_up)
    
    deallocate( spectra_up, spectra_dn )
    
    call out_data_sub(opt//'-dn.dat', data_dn)
    call out_data_sub(opt//'-up.dat', data_up)

    open(unit=8, file='ice_ranges_'//opt//'_dn_max', status='new', action='write')
      write(8,*) ceiling( maxval(data_dn) )
    close(8)
    open(unit=8, file='ice_ranges_'//opt//'_dn_min', status='new', action='write')
      write(8,*) floor( minval(data_dn) )
    close(8)
    open(unit=8, file='ice_ranges_'//opt//'_dn_cnt', status='new', action='write')
      write(8,*) ceiling( ( (( maxval(data_dn) - minval(data_dn) ) / 8 ) / 50 ) ) * 50
    close(8)

    open(unit=8, file='ice_ranges_'//opt//'_up_max', status='new', action='write')
      write(8,*) ceiling( maxval(data_up) )
    close(8)
    open(unit=8, file='ice_ranges_'//opt//'_up_min', status='new', action='write')
      write(8,*) floor( minval(data_up) )
    close(8)
    open(unit=8, file='ice_ranges_'//opt//'_up_cnt', status='new', action='write')
      write(8,*) ceiling( ( ( ( maxval(data_up) - minval(data_up) ) / 8 ) / 50 ) ) * 50
    close(8)
    
    deallocate( data_up, data_dn )
    
  end subroutine harm_analysis_ice_sub

    subroutine get_spectra_sub(path)
      character(len=*), intent(in) :: path
      character(len=10)            :: subor
      integer                      :: j, m, n, error
      complex(kind=dbl)            :: u_dn, u_up
    
      !Najdi najnovsi subor s tvarovymi koeficientami
      n = 0
      do
        write(subor, '(1I4)') n
    
        open(unit=8, file=path//trim(adjustl(subor))//'.dat', status='old', action='read', iostat=error)
        if (error /= 0) then
          write(subor, '(1I4)') n-1
          exit
        end if
        close(8)
    
        n = n + 1
      end do
      
      !Nacitaj spektrum najnovsich koeficientov
      open(unit=8, file=path//trim(adjustl(subor))//'.dat', status='old', action='read')
    
      do
        read(8,*,iostat=error) j, m, u_dn, u_up; if (error /= 0) exit
        
        spectra_dn(j*(j+1)/2+m+1) = u_dn
        spectra_up(j*(j+1)/2+m+1) = u_up
      end do
    
      close(8)
      
    end subroutine get_spectra_sub
    
end module OutputIceMod
