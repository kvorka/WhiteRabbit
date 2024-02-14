module SphericalHarmonics
  use Clebsch_Gordan
  use Legendre_polynomials
  use Fourier_transform
  use sph_vectors
  implicit none
  
  type, public :: T_lateralGrid
    integer,                     private :: jmax, jms, jms1, jms2, jmv, jmv1, maxj, nLegendre, nFourier
    real(kind=dbl),              private :: tolm, scale
    real(kind=dbl), allocatable, private :: roots(:), fftLege(:), amjrr(:), bmjrr(:)
    type(T_fft),                 private :: fourtrans
    
    contains
    
    procedure :: init_sub       => init_harmonics_sub
    procedure :: deallocate_sub => deallocate_harmonics_sub
    procedure :: partial_backward_2_sub, partial_forward_2_sub
    procedure :: partial_backward_4_sub, partial_forward_4_sub
    procedure :: partial_backward_8_sub, partial_forward_8_sub
    procedure :: partial_backward_16_sub, partial_forward_16_sub
    procedure :: grid_op_2_vcsum_sub, grid_op_4_vcsum_sub, grid_op_8_vcsum_sub, grid_op_16_vcsum_sub
    procedure :: grid_op_2_vcvv_sub, grid_op_4_vcvv_sub, grid_op_8_vcvv_sub, grid_op_16_vcvv_sub
    procedure :: grid_op_2_vcvgv_sub, grid_op_4_vcvgv_sub, grid_op_8_vcvgv_sub, grid_op_16_vcvgv_sub
    procedure :: grid_op_2_vcvv_vcvgv_sub, grid_op_4_vcvv_vcvgv_sub, grid_op_8_vcvv_vcvgv_sub, grid_op_16_vcvv_vcvgv_sub
    procedure :: vctol_sub, vcsum_sub, vcvv_sub, vcvgv_sub, vcvv_vcvgv_sub
    
  end type T_lateralGrid
  
  interface
    module subroutine init_harmonics_sub(this, jmax)
      class(T_lateralGrid), intent(inout) :: this
      integer,              intent(in)    :: jmax
    end subroutine init_harmonics_sub
    
    module pure subroutine deallocate_harmonics_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_harmonics_sub
    
    module pure subroutine partial_backward_2_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: n
      real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
      real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    end subroutine partial_backward_2_sub
    
    module pure subroutine partial_forward_2_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: n
      real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
      real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
      complex(kind=dbl),    intent(inout) :: cr(*)
      complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    end subroutine partial_forward_2_sub
    
    module pure subroutine partial_backward_4_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: n
      real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
      real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    end subroutine partial_backward_4_sub
    
    module pure subroutine partial_forward_4_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: n
      real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
      real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
      complex(kind=dbl),    intent(inout) :: cr(*)
      complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    end subroutine partial_forward_4_sub
    
    module pure subroutine partial_backward_8_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: n
      real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
      real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    end subroutine partial_backward_8_sub
    
    module pure subroutine partial_forward_8_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: n
      real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
      real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
      complex(kind=dbl),    intent(inout) :: cr(*)
      complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    end subroutine partial_forward_8_sub
    
    module pure subroutine partial_backward_16_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: n
      real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
      real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    end subroutine partial_backward_16_sub
    
    module pure subroutine partial_forward_16_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: n
      real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
      real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
      complex(kind=dbl),    intent(inout) :: cr(*)
      complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    end subroutine partial_forward_16_sub
    
    module subroutine vctol_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine vctol_sub
    
    module pure subroutine vcsum_sub(this, cajm, cbjm, cjm)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcsum_sub
    
    module pure subroutine vcvv_sub(this, cajml, cbjml, cjm)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cajml(*), cbjml(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcvv_sub
    
    module pure subroutine vcvgv_sub(this, ri, dv_r, v, cjm)
      class(T_lateralGrid), intent(in)  :: this
      real(kind=dbl),       intent(in)  :: ri
      complex(kind=dbl),    intent(in)  :: v(*), dv_r(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcvgv_sub
    
    module pure subroutine vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
      class(T_lateralGrid), intent(in)  :: this
      real(kind=dbl),       intent(in)  :: ri
      complex(kind=dbl),    intent(in)  :: dv_r(*), q(*), v(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcvv_vcvgv_sub
    
    module pure subroutine grid_op_2_vcsum_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_2_vcsum_sub
    
    module pure subroutine grid_op_4_vcsum_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_4_vcsum_sub
    
    module pure subroutine grid_op_8_vcsum_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_8_vcsum_sub
    
    module pure subroutine grid_op_16_vcsum_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_16_vcsum_sub
    
    module pure subroutine grid_op_2_vcvv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_2_vcvv_sub
    
    module pure subroutine grid_op_4_vcvv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_4_vcvv_sub
    
    module pure subroutine grid_op_8_vcvv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_8_vcvv_sub
    
    module pure subroutine grid_op_16_vcvv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_16_vcvv_sub
    
    module pure subroutine grid_op_2_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_2_vcvgv_sub
    
    module pure subroutine grid_op_4_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_4_vcvgv_sub
    
    module pure subroutine grid_op_8_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_8_vcvgv_sub
    
    module pure subroutine grid_op_16_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_16_vcvgv_sub
    
    module pure subroutine grid_op_2_vcvv_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_2_vcvv_vcvgv_sub
    
    module pure subroutine grid_op_4_vcvv_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_4_vcvv_vcvgv_sub
    
    module pure subroutine grid_op_8_vcvv_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_8_vcvv_vcvgv_sub
    
    module pure subroutine grid_op_16_vcvv_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(out)   :: grid(*)
      complex(kind=dbl),    intent(inout) :: sumNS(*)
    end subroutine grid_op_16_vcvv_vcvgv_sub
  end interface
  
end module SphericalHarmonics
