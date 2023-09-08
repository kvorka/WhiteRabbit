module SphericalHarmonics
  use Clebsch_Legendre
  use, intrinsic :: iso_c_binding
  implicit none
  include '/usr/include/fftw3.f03'

  integer, parameter :: fftw_flags = FFTW_EXHAUSTIVE + FFTW_UNALIGNED + FFTW_DESTROY_INPUT
  integer, parameter :: step       = 8
  
  type, public :: T_lateralGrid
    integer,                     private :: jmax, jms, jms1, jms2, jmv, jmv1, maxj, nLegendre, nFourier
    real(kind=dbl), allocatable, private :: roots(:), fftLege(:), ish(:), amjrr(:), bmjrr(:)
    type(C_ptr),                 private :: fftw_06_c2r, fftw_01_r2c
    type(C_ptr),                 private :: fftw_16_c2r, fftw_03_r2c
    type(C_ptr),                 private :: fftw_19_c2r, fftw_04_r2c
  
    contains
  
    procedure :: init_sub       => init_harmonics_sub
    procedure :: deallocate_sub => deallocate_harmonics_sub

    procedure :: init_vcvv_sub
    procedure :: vcvv_fn
    procedure :: deallocate_fftw_vcvv_sub

    procedure :: init_vcsv_vcvgv_sub
    procedure :: vcsv_vcvgv_fn
    procedure :: deallocate_fftw_vcsv_vcvgv_sub

    procedure :: init_vcsv_vcvv_vcvgv_sub
    procedure :: vcsv_vcvv_vcvgv_sub
    procedure :: deallocate_fftw_vcsv_vcvv_vcvgv_sub

  end type T_lateralGrid
  
  interface
    module subroutine init_harmonics_sub(this, jmax)
      class(T_lateralGrid), intent(inout) :: this
      integer,              intent(in)    :: jmax
    end subroutine init_harmonics_sub
    
    module subroutine deallocate_harmonics_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_harmonics_sub
    
    module subroutine init_vcvv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine init_vcvv_sub
    
    module function vcvv_fn(this, cajml, cbjml) result(cjm)
      class(T_lateralGrid), intent(in) :: this
      complex(kind=dbl),    intent(in) :: cajml(:), cbjml(:)
      complex(kind=dbl)                :: cjm(this%jms)
    end function vcvv_fn
    
    module subroutine deallocate_fftw_vcvv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_fftw_vcvv_sub
    
    module subroutine init_vcsv_vcvgv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine init_vcsv_vcvgv_sub
    
    module function vcsv_vcvgv_fn(this, ri, dv_r, v) result(cjml)
      class(T_lateralGrid), intent(in) :: this
      real(kind=dbl),       intent(in) :: ri
      complex(kind=dbl),    intent(in) :: v(:), dv_r(:)
      complex(kind=dbl)                :: cjml(this%jmv)
    end function vcsv_vcvgv_fn
    
    module subroutine deallocate_fftw_vcsv_vcvgv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_fftw_vcsv_vcvgv_sub
    
    module subroutine init_vcsv_vcvv_vcvgv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine init_vcsv_vcvv_vcvgv_sub
    
    module subroutine vcsv_vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
      class(T_lateralGrid), intent(in)  :: this
      real(kind=dbl),       intent(in)  :: ri
      complex(kind=dbl),    intent(in)  :: dv_r(:), q(:), v(:)
      complex(kind=dbl),    intent(out) :: cjm(:,:)
    end subroutine vcsv_vcvv_vcvgv_sub
    
    module subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub
  end interface
  
end module SphericalHarmonics
