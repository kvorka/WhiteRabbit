module OutputIceMod
  use IceConstants
  use OutputMod
  implicit none
  
  public :: harm_analysis_ice_sub
  public :: zonal_analysis_ice_sub
  public :: save_spectra_ice_sub

  private :: cnt_fn
 
  contains

  subroutine zonal_analysis_ice_sub(opt)
    character(len=*),   intent(in) :: opt
    integer                        :: jm, jms, error, j, m
    complex(kind=dbl), allocatable :: spectra_up(:), spectra_dn(:)
    real(kind=dbl),    allocatable :: data_up(:,:), data_dn(:,:)
    
    jms = jmax_ice*(jmax_ice+1)/2+jmax_ice+1; allocate( spectra_up(jms), spectra_dn(jms) )
    
    open(unit=1, file=opt//'-dn.spec', status='old', action='read')
    open(unit=2, file=opt//'-up.spec', status='old', action='read')
      do
        read(1,*,iostat=error) jm, spectra_dn(jm) ; if (error /= 0) exit
        read(2,*) jm, spectra_up(jm)
      end do
    close(1)
    close(2)

    do j = 0, jmax_ice
      do m = 1, j
        jm = j*(j+1)/2+m+1

        spectra_dn(jm) = cmplx(0._dbl, 0._dbl, kind=dbl)
        spectra_up(jm) = cmplx(0._dbl, 0._dbl, kind=dbl)
      end do
    end do
    
    allocate( data_dn(2*nth,nth) ); data_dn = 0._dbl; call toGrid_sub(jmax_ice, spectra_dn, data_dn); deallocate(spectra_dn)
      call out_zondata_sub(opt//'-zon-dn.dat', data_dn)
    deallocate(data_dn)
    
    allocate( data_up(2*nth,nth) ); data_up = 0._dbl; call toGrid_sub(jmax_ice, spectra_up, data_up); deallocate(spectra_up)
      call out_zondata_sub(opt//'-zon-up.dat', data_up)
    deallocate(data_up)

  end subroutine zonal_analysis_ice_sub

  subroutine harm_analysis_ice_sub(opt)
    character(len=*),   intent(in) :: opt
    integer                        :: jm, jms, error
    complex(kind=dbl), allocatable :: spectra_up(:), spectra_dn(:)
    real(kind=dbl),    allocatable :: data_up(:,:), data_dn(:,:)
    
    jms = jmax_ice*(jmax_ice+1)/2+jmax_ice+1; allocate( spectra_up(jms), spectra_dn(jms) )
    
    open(unit=1, file=opt//'-dn.spec', status='old', action='read')
    open(unit=2, file=opt//'-up.spec', status='old', action='read')
      do
        read(1,*,iostat=error) jm, spectra_dn(jm) ; if (error /= 0) exit
        read(2,*) jm, spectra_up(jm)
      end do
    close(1)
    close(2)
    
    allocate( data_dn(2*nth,nth) ); data_dn = 0._dbl; call toGrid_sub(jmax_ice, spectra_dn, data_dn); deallocate(spectra_dn)
      call out_data_sub(opt//'-dn.dat', data_dn)
        open(unit=8, file='i_'//opt//'_dn_max', status='new', action='write'); write(8,*) ceiling(maxval(data_dn)); close(8)
        open(unit=8, file='i_'//opt//'_dn_min', status='new', action='write'); write(8,*) floor(minval(data_dn))  ; close(8)
        open(unit=8, file='i_'//opt//'_dn_cnt', status='new', action='write'); write(8,*) cnt_fn(data_dn)         ; close(8)
    deallocate(data_dn)
    
    allocate( data_up(2*nth,nth) ); data_up = 0._dbl; call toGrid_sub(jmax_ice, spectra_up, data_up); deallocate(spectra_up)
      call out_data_sub(opt//'-up.dat', data_up)
        open(unit=8, file='i_'//opt//'_up_max', status='new', action='write'); write(8,*) ceiling(maxval(data_up)); close(8)
        open(unit=8, file='i_'//opt//'_up_min', status='new', action='write'); write(8,*) floor(minval(data_up))  ; close(8)
        open(unit=8, file='i_'//opt//'_up_cnt', status='new', action='write'); write(8,*) cnt_fn(data_up)         ; close(8)
    deallocate(data_up)
    
  end subroutine harm_analysis_ice_sub

    integer function cnt_fn(data_in)
      real(kind=dbl), intent(in) :: data_in(:,:)
      
      cnt_fn = ceiling( ( (( maxval(data_in) - minval(data_in) ) / 8 ) / 50 ) ) * 50
      
    end function cnt_fn

  subroutine save_spectra_ice_sub(path, opt)
    character(len=*),   intent(in) :: path, opt
    integer                        :: j, m, jms, n, error
    complex(kind=dbl)              :: u_dn, u_up
    complex(kind=dbl), allocatable :: spectra_up(:), spectra_dn(:)
    
    jms = jmax_ice*(jmax_ice+1)/2+jmax_ice+1 ; allocate( spectra_up(jms), spectra_dn(jms) )
    
    n = 0
    do
      open(unit=8, file=path//trim(adjustl(int2str_fn(n)))//'.dat', status='old', action='read', iostat=error)
        if (error /= 0) then
          n = n - 1
          exit
        end if
      close(8)
      
      n = n + 1
    end do
    
    open(unit=8, file=path//trim(adjustl(int2str_fn(n)))//'.dat', status='old', action='read')
      do
        read(8,*,iostat=error) j, m, u_dn, u_up; if (error /= 0) exit
        
        spectra_dn(j*(j+1)/2+m+1) = u_dn
        spectra_up(j*(j+1)/2+m+1) = u_up
      end do
    close(8)
    
    call out_spectra1_sub(opt//'-dn.spec', spectra_dn)
    call out_spectra1_sub(opt//'-up.spec', spectra_up)
    
    deallocate( spectra_up, spectra_dn )
    
  end subroutine save_spectra_ice_sub
    
end module OutputIceMod
