module SphericalHarmonics
  use Clebsch_Legendre
  use FFT_mod
  implicit none
  
  integer, parameter :: step = 8
  
  type, public :: T_lateralGrid
    integer,                     private :: jmax, jms, jms1, jms2, jmv, jmv1, maxj, nLegendre, nFourier
    real(kind=dbl),              private :: tolm
    real(kind=dbl), allocatable, private :: roots(:), fftLege(:), amjrr(:), bmjrr(:)
    type(T_fft),                 private :: fourtrans
    
    contains
    
    procedure :: init_sub       => init_harmonics_sub
    procedure :: deallocate_sub => deallocate_harmonics_sub
    procedure :: vcvv_sub, vcsv_vcvgv_sub, vcsv_vcvv_vcvgv_sub
    
  end type T_lateralGrid
  
  interface
    module subroutine init_harmonics_sub(this, jmax)
      class(T_lateralGrid), intent(inout) :: this
      integer,              intent(in)    :: jmax
    end subroutine init_harmonics_sub
    
    module subroutine deallocate_harmonics_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_harmonics_sub
    
    module subroutine vcvv_sub(this, cajml, cbjml, cjm)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cajml(:), cbjml(:)
      complex(kind=dbl),    intent(out) :: cjm(:)
    end subroutine vcvv_sub
    
    module subroutine vcsv_vcvgv_sub(this, ri, dv_r, v, cjm)
      class(T_lateralGrid), intent(in)  :: this
      real(kind=dbl),       intent(in)  :: ri
      complex(kind=dbl),    intent(in)  :: v(:), dv_r(:)
      complex(kind=dbl),    intent(out) :: cjm(:,:)
    end subroutine vcsv_vcvgv_sub
    
    module pure subroutine vcsv_vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
      class(T_lateralGrid), intent(in)  :: this
      real(kind=dbl),       intent(in)  :: ri
      complex(kind=dbl),    intent(in)  :: dv_r(:), q(:), v(:)
      complex(kind=dbl),    intent(out) :: cjm(:,:)
    end subroutine vcsv_vcvv_vcvgv_sub
  end interface
  
end module SphericalHarmonics
