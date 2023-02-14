module OutputOceanMod
  use OceanConstants
  use OutputMod
  use Vector_analysis
  implicit none
  
  integer,        parameter :: n_out = 120
  real(kind=dbl), parameter :: treshold = 0.8_dbl
  integer,        parameter :: i1    = 200
  integer,        parameter :: i2    = 400

  integer, parameter :: jms  = jmax_ocean * (jmax_ocean+1) / 2 + jmax_ocean + 1
  integer, parameter :: jms1 = (jmax_ocean+1)*(jmax_ocean+2)/2 + (jmax_ocean+1) + 1
  integer, parameter :: jmv  = 3 * (jmax_ocean*(jmax_ocean+1)/2 + jmax_ocean) + 1
  
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
    module subroutine harm_analysis_temp_sub(path)
      character(len=*), intent(in) :: path
    end subroutine harm_analysis_temp_sub

    module subroutine save_spectra_temp_sub(path)
      character(len=*), intent(in) :: path
    end subroutine save_spectra_temp_sub
  end interface

  interface
    module subroutine harm_analysis_flux_sub(path)
      character(len=*), intent(in) :: path
    end subroutine harm_analysis_flux_sub

    module subroutine save_spectra_flux_sub(path)
      character(len=*), intent(in) :: path
    end subroutine save_spectra_flux_sub
  end interface

  interface
    module subroutine harm_analysis_rad_velc_sub(path)
      character(len=*), intent(in) :: path
    end subroutine harm_analysis_rad_velc_sub

    module subroutine harm_analysis_zon_velc_sub(path)
      character(len=*), intent(in) :: path
    end subroutine harm_analysis_zon_velc_sub

    module subroutine save_spectra_velc_sub(path)
      character(len=*), intent(in) :: path
    end subroutine save_spectra_velc_sub
  end interface

  interface
    module subroutine nuss_curve_sub(path)
      character(len=*), intent(in) :: path
    end subroutine nuss_curve_sub
  end interface

  interface
    module subroutine get_zonal_sub(data_in, data_out)
      real(kind=dbl), intent(in)  :: data_in(:,:)
      real(kind=dbl), intent(out) :: data_out(:)
    end subroutine get_zonal_sub

    module subroutine load_data_sub(dim_in, file_in, data_out)
      character(len=*),  intent(in)  :: file_in
      integer,           intent(in)  :: dim_in
      complex(kind=dbl), intent(out) :: data_out(:,:)
    end subroutine load_data_sub

    module subroutine save_data_sub(file_in, data_in)
      character(len=*), intent(in) :: file_in
      real(kind=dbl),   intent(in) :: data_in(:,:)
    end subroutine save_data_sub
  end interface

end module OutputOceanMod
  
