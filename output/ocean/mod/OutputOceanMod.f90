module OutputOceanMod
  use OceanConstants
  use Paths
  use OutputMod
  use Vector_analysis
  use Spherical_func
  implicit none
  
  integer, parameter, public :: n_out = 120
  integer, parameter, public :: jms  =    (jmax_ocean  )*(jmax_ocean+1)/2 + (jmax_ocean  )  + 1
  integer, parameter, public :: jms1 =    (jmax_ocean+1)*(jmax_ocean+2)/2 + (jmax_ocean+1)  + 1
  integer, parameter, public :: jmv  = 3*((jmax_ocean  )*(jmax_ocean+1)/2 + (jmax_ocean  )) + 1
  
  public :: nuss_curve_sub
  
  public :: harm_analysis_flux_sub
  public :: harm_analysis_temp_sub
  public :: harm_analysis_rad_velc_sub
  public :: harm_analysis_zon_velc_sub
  
  public :: save_spectra_velc_sub
  public :: save_spectra_temp_sub
  public :: save_spectra_flux_sub
  
  public :: get_zonal_sub
  public :: load_data_sub
  public :: save_data_sub
  
  interface
    module subroutine nuss_curve_sub()
    end subroutine nuss_curve_sub
    
    module subroutine harm_analysis_flux_sub()
    end subroutine harm_analysis_flux_sub
    
    module subroutine save_spectra_flux_sub()
    end subroutine save_spectra_flux_sub
    
    module subroutine harm_analysis_temp_sub()
    end subroutine harm_analysis_temp_sub
    
    module subroutine save_spectra_temp_sub()
    end subroutine save_spectra_temp_sub
    
    module subroutine harm_analysis_rad_velc_sub()
    end subroutine harm_analysis_rad_velc_sub
    
    module subroutine harm_analysis_zon_velc_sub()
    end subroutine harm_analysis_zon_velc_sub
    
    module subroutine save_spectra_velc_sub()
    end subroutine save_spectra_velc_sub
  end interface
  
  contains
  
  subroutine get_zonal_sub(data_in, data_out)
    real(kind=dbl), intent(in)  :: data_in(:,:)
    real(kind=dbl), intent(out) :: data_out(:)
    integer                     :: ith
    
    data_out(1) = data_in(1,1)
    
    do ith = 1, nth-1
      data_out(ith+1) = ( data_in(1,ith) + data_in(1,ith+1) ) / 2
    end do
    
    data_out(nth+1) = data_in(1,nth)
    
  end subroutine get_zonal_sub
  
  subroutine load_data_sub(dim_in, file_in, data_out)
    character(len=*),  intent(in)  :: file_in
    integer,           intent(in)  :: dim_in
    complex(kind=dbl), intent(out) :: data_out(:,:)
    real(kind=dbl),    allocatable :: r_in(:)
    complex(kind=dbl), allocatable :: data_in(:,:)
    integer                        :: i, ii
    real(kind=dbl)                 :: r

    allocate( r_in(nd_ocean+1), data_in(dim_in, nd_ocean+1) )

      open(unit=1, file=file_in, status='old', action='read')
        do i = 1, nd_ocean+1
          read(1,*) r_in(i), data_in(:,i)
        end do
      close(1)

      do i = 1, n_out
        r = r_ud_ocean/(1-r_ud_ocean) + (i-1._dbl)/(n_out-1._dbl)

        do ii = 1, nd_ocean
          if ( (r >= r_in(ii)) .and. (r <= r_in(ii+1)) ) then
            data_out(:,i) = ( (r-r_in(ii)) * data_in(:,ii+1) + (r_in(ii+1)-r) * data_in(:,ii) ) / ( r_in(ii+1)-r_in(ii) )
            exit
          end if
        end do
      end do
    
    deallocate( r_in, data_in )
    
    data_out(1,:) = czero

  end subroutine load_data_sub

  subroutine save_data_sub(file_in, data_in)
    character(len=*), intent(in) :: file_in
    real(kind=dbl),   intent(in) :: data_in(:,:)
    integer                      :: i, k

    open(unit=2, file=file_in, status='new', action='write')
    
      do k = 0, nth
        do i = 1, n_out
          write(2,*) k+90, 480+i, data_in(k+1,i)
        end do
      end do

    close(2)

  end subroutine save_data_sub

end module OutputOceanMod
  
