module SphericalHarmonics
  use Nulify
  use Clebsch_Gordan
  use Legendre_function
  use Fourier_transform
  implicit none
  
  integer, parameter :: addmissible_jmax(47) = [   5,   7,   9,  13,  15,  21,  27,  29,  33,  37,  45, 47,  51,  57,  61,  69,  &
                                               &  77,  87,  93,  97, 105, 117, 125, 141, 147, 157, 159, 177, 189, 197, 213, 237, &
                                               & 247, 253, 267, 285, 297, 317, 321, 357, 381, 397, 429, 447, 477, 483, 497       ]
  
  type, public :: T_lateralGrid
    integer,                     private :: jmax, jmax1, jmax2, jmax3
    integer,                     private :: jms, jms1, jms2
    integer,                     private :: jmv, jmv1
    integer,                     private :: jmt
    integer,                     private :: nLegendre, nFourier
    integer,        allocatable, private :: maxm(:)
    real(kind=dbl),              private :: tolm, scale
    real(kind=dbl), allocatable, private :: roots(:), fftLege(:), amjrr(:), bmjrr(:), pmm(:,:)
    type(T_fft),                 private :: fourtrans
    
    contains
    
    procedure, pass   :: init_sub       => init_harmonics_sub
    procedure, pass   :: deallocate_sub => deallocate_harmonics_sub
    procedure, pass   :: scal2scal_jm_to_mj_sub, scal2scal_mj_to_jm_sub
    procedure, pass   :: vec2vec_jml_to_jml_sub, vec2scal_jml_to_mj_sub, scal2vecscal_mj_to_jm_sub
    procedure, pass   :: gradvec2vec_jmlk_to_jml_sub, devtens2scal_jml2_to_mj_sub, scal2devtens_mj_to_jml2_sub
    procedure, pass   :: rescale_sub, get_maxm_fn, vctol_sub
    procedure, pass   :: pmj_backward_set_2_sub, pmj_backward_set_4_sub, pmj_backward_set_8_sub, pmj_backward_set_16_sub
    procedure, pass   :: pmj_backward_rec_2_sub, pmj_backward_rec_4_sub, pmj_backward_rec_8_sub, pmj_backward_rec_16_sub
    procedure, pass   :: pmj_forward_set_2_sub, pmj_forward_set_4_sub, pmj_forward_set_8_sub, pmj_forward_set_16_sub
    procedure, pass   :: pmj_forward_rec_2_sub, pmj_forward_rec_4_sub, pmj_forward_rec_8_sub, pmj_forward_rec_16_sub
    procedure, nopass :: pmj_backward_recomb_2_sub, pmj_backward_recomb_4_sub, pmj_backward_recomb_8_sub, pmj_backward_recomb_16_sub
    procedure, nopass :: pmj_forward_recomb_2_sub, pmj_forward_recomb_4_sub, pmj_forward_recomb_8_sub, pmj_forward_recomb_16_sub
    procedure, pass   :: lege_transform_sub
    procedure, pass   :: grid_op_2_vcsum_sub, grid_op_4_vcsum_sub, grid_op_8_vcsum_sub, grid_op_16_vcsum_sub
    procedure, pass   :: grid_op_2_vcst_sub, grid_op_4_vcst_sub, grid_op_8_vcst_sub, grid_op_16_vcst_sub
    procedure, pass   :: grid_op_2_vcvv_sub, grid_op_4_vcvv_sub, grid_op_8_vcvv_sub, grid_op_16_vcvv_sub
    procedure, pass   :: grid_op_2_vcvgv_sub, grid_op_4_vcvgv_sub, grid_op_8_vcvgv_sub, grid_op_16_vcvgv_sub
    procedure, pass   :: grid_op_2_vcvv_vcvgv_sub, grid_op_4_vcvv_vcvgv_sub, grid_op_8_vcvv_vcvgv_sub, grid_op_16_vcvv_vcvgv_sub
    procedure, pass   :: vcsum_sub, vcst_sub, vcvv_sub, vcvgv_sub, vcvv_vcvgv_sub
    
  end type T_lateralGrid
  
  interface
    module subroutine init_harmonics_sub(this, jmax)
      class(T_lateralGrid), intent(inout) :: this
      integer,              intent(in)    :: jmax
    end subroutine init_harmonics_sub
    
    module pure subroutine deallocate_harmonics_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_harmonics_sub
    
    module pure subroutine rescale_sub( this, arr, length )
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: length
      complex(kind=dbl),    intent(inout) :: arr(*)
    end subroutine rescale_sub
    
    module pure integer function get_maxm_fn(this, i, i2)
      class(T_lateralGrid), intent(in) :: this
      integer,              intent(in) :: i, i2
    end function get_maxm_fn
    
    module pure subroutine pmj_backward_set_2_sub(this, i, m, nsum, pmj2, pmj1, pmj, cc, legesum)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: nsum, i, m
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: legesum(2,*)
    end subroutine pmj_backward_set_2_sub
    
    module pure subroutine pmj_forward_set_2_sub(this, i, m, nsum, pmj2, pmj1, pmj, legesum, cr)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: i, m, nsum
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(inout) :: legesum(2,*)
      complex(kind=dbl),    intent(inout) :: cr(*)
    end subroutine pmj_forward_set_2_sub
    
    module pure subroutine pmj_backward_set_4_sub(this, i, m, nsum, pmj2, pmj1, pmj, cc, legesum)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: i, m, nsum
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: legesum(4,*)
    end subroutine pmj_backward_set_4_sub
    
    module pure subroutine pmj_forward_set_4_sub(this, i, m, nsum, pmj2, pmj1, pmj, legesum, cr)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: i, m, nsum
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(inout) :: legesum(4,*)
      complex(kind=dbl),    intent(inout) :: cr(*)
    end subroutine pmj_forward_set_4_sub
    
    module pure subroutine pmj_backward_set_8_sub(this, i, m, nsum, pmj2, pmj1, pmj, cc, legesum)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: i, m, nsum
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: legesum(8,*)
    end subroutine pmj_backward_set_8_sub
    
    module pure subroutine pmj_forward_set_8_sub(this, i, m, nsum, pmj2, pmj1, pmj, legesum, cr)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: i, m, nsum
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(inout) :: legesum(8,*)
      complex(kind=dbl),    intent(inout) :: cr(*)
    end subroutine pmj_forward_set_8_sub
    
    module pure subroutine pmj_backward_set_16_sub(this, i, m, nsum, pmj2, pmj1, pmj, cc, legesum)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: i, m, nsum
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: legesum(16,*)
    end subroutine pmj_backward_set_16_sub
    
    module pure subroutine pmj_forward_set_16_sub(this, i, m, nsum, pmj2, pmj1, pmj, legesum, cr)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: i, m, nsum
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(inout) :: legesum(16,*)
      complex(kind=dbl),    intent(inout) :: cr(*)
    end subroutine pmj_forward_set_16_sub
    
    module pure subroutine pmj_backward_rec_2_sub(this, mj, nsum, cosx, pmj2, pmj1, pmj, cc, legesum)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: mj, nsum
      real(kind=dbl),       intent(in)    :: cosx(*)
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: legesum(2,*)
    end subroutine pmj_backward_rec_2_sub
    
    module pure subroutine pmj_forward_rec_2_sub(this, mj, nsum, cosx, pmj2, pmj1, pmj, legesum, cr)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: mj, nsum
      real(kind=dbl),       intent(in)    :: cosx(*)
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: legesum(2,*)
      complex(kind=dbl),    intent(inout) :: cr(*)
    end subroutine pmj_forward_rec_2_sub
    
    module pure subroutine pmj_backward_rec_4_sub(this, mj, nsum, cosx, pmj2, pmj1, pmj, cc, legesum)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: mj, nsum
      real(kind=dbl),       intent(in)    :: cosx(*)
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: legesum(4,*)
    end subroutine pmj_backward_rec_4_sub
    
    module pure subroutine pmj_forward_rec_4_sub(this, mj, nsum, cosx, pmj2, pmj1, pmj, legesum, cr)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: mj, nsum
      real(kind=dbl),       intent(in)    :: cosx(*)
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: legesum(4,*)
      complex(kind=dbl),    intent(inout) :: cr(*)
    end subroutine pmj_forward_rec_4_sub
    
    module pure subroutine pmj_backward_rec_8_sub(this, mj, nsum, cosx, pmj2, pmj1, pmj, cc, legesum)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: mj, nsum
      real(kind=dbl),       intent(in)    :: cosx(*)
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: legesum(8,*)
    end subroutine pmj_backward_rec_8_sub
    
    module pure subroutine pmj_forward_rec_8_sub(this, mj, nsum, cosx, pmj2, pmj1, pmj, legesum, cr)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: mj, nsum
      real(kind=dbl),       intent(in)    :: cosx(*)
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: legesum(8,*)
      complex(kind=dbl),    intent(inout) :: cr(*)
    end subroutine pmj_forward_rec_8_sub
    
    module pure subroutine pmj_backward_rec_16_sub(this, mj, nsum, cosx, pmj2, pmj1, pmj, cc, legesum)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: mj, nsum
      real(kind=dbl),       intent(in)    :: cosx(*)
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: cc(*)
      complex(kind=dbl),    intent(inout) :: legesum(16,*)
    end subroutine pmj_backward_rec_16_sub
    
    module pure subroutine pmj_forward_rec_16_sub(this, mj, nsum, cosx, pmj2, pmj1, pmj, legesum, cr)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: mj, nsum
      real(kind=dbl),       intent(in)    :: cosx(*)
      real(kind=dbl),       intent(inout) :: pmj1(*), pmj2(*), pmj(*)
      complex(kind=dbl),    intent(in)    :: legesum(16,*)
      complex(kind=dbl),    intent(inout) :: cr(*)
    end subroutine pmj_forward_rec_16_sub
    
    module pure subroutine pmj_backward_recomb_2_sub(nsum, ssym, asym, sumN, sumS)
      integer,              intent(in)    :: nsum
      complex(kind=dbl),    intent(in)    :: ssym(2,*), asym(2,*)
      complex(kind=dbl),    intent(inout) :: sumN(nsum,*), sumS(nsum,*)
    end subroutine pmj_backward_recomb_2_sub
    
    module pure subroutine pmj_forward_recomb_2_sub(nsum, weight, sumN, sumS, ssym, asym)
      integer,           intent(in)  :: nsum
      real(kind=dbl),    intent(in)  :: weight(*)
      complex(kind=dbl), intent(in)  :: sumN(nsum,*), sumS(nsum,*)
      complex(kind=dbl), intent(out) :: ssym(2,*), asym(2,*)
    end subroutine pmj_forward_recomb_2_sub
    
    module pure subroutine pmj_backward_recomb_4_sub(nsum, ssym, asym, sumN, sumS)
      integer,              intent(in)    :: nsum
      complex(kind=dbl),    intent(in)    :: ssym(4,*), asym(4,*)
      complex(kind=dbl),    intent(inout) :: sumN(nsum,*), sumS(nsum,*)
    end subroutine pmj_backward_recomb_4_sub
    
    module pure subroutine pmj_forward_recomb_4_sub(nsum, weight, sumN, sumS, ssym, asym)
      integer,           intent(in)  :: nsum
      real(kind=dbl),    intent(in)  :: weight(*)
      complex(kind=dbl), intent(in)  :: sumN(nsum,*), sumS(nsum,*)
      complex(kind=dbl), intent(out) :: ssym(4,*), asym(4,*)
    end subroutine pmj_forward_recomb_4_sub
    
    module pure subroutine pmj_backward_recomb_8_sub(nsum, ssym, asym, sumN, sumS)
      integer,              intent(in)    :: nsum
      complex(kind=dbl),    intent(in)    :: ssym(8,*), asym(8,*)
      complex(kind=dbl),    intent(inout) :: sumN(nsum,*), sumS(nsum,*)
    end subroutine pmj_backward_recomb_8_sub
    
    module pure subroutine pmj_forward_recomb_8_sub(nsum, weight, sumN, sumS, ssym, asym)
      integer,           intent(in)  :: nsum
      real(kind=dbl),    intent(in)  :: weight(*)
      complex(kind=dbl), intent(in)  :: sumN(nsum,*), sumS(nsum,*)
      complex(kind=dbl), intent(out) :: ssym(8,*), asym(8,*)
    end subroutine pmj_forward_recomb_8_sub
    
    module pure subroutine pmj_backward_recomb_16_sub(nsum, ssym, asym, sumN, sumS)
      integer,              intent(in)    :: nsum
      complex(kind=dbl),    intent(in)    :: ssym(16,*), asym(16,*)
      complex(kind=dbl),    intent(inout) :: sumN(nsum,*), sumS(nsum,*)
    end subroutine pmj_backward_recomb_16_sub
    
    module pure subroutine pmj_forward_recomb_16_sub(nsum, weight, sumN, sumS, ssym, asym)
      integer,           intent(in)  :: nsum
      real(kind=dbl),    intent(in)  :: weight(*)
      complex(kind=dbl), intent(in)  :: sumN(nsum,*), sumS(nsum,*)
      complex(kind=dbl), intent(out) :: ssym(16,*), asym(16,*)
    end subroutine pmj_forward_recomb_16_sub
    
    module pure subroutine lege_transform_sub(this, nforw, nback, cc, cr, grid_2_sub, grid_4_sub, grid_8_sub, grid_16_sub)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: nforw, nback
      complex(kind=dbl),    intent(in)    :: cc(nback,*)
      complex(kind=dbl),    intent(inout) :: cr(nforw,*)
      
      interface
        pure subroutine grid_2_sub(sph, gxyz, sumNS)
          import dbl, T_lateralGrid
          class(T_lateralGrid),   intent(in)    :: sph
          real(kind=dbl), target, intent(out)   :: gxyz(*)
          complex(kind=dbl),      intent(inout) :: sumNS(*)
        end subroutine grid_2_sub
        
        pure subroutine grid_4_sub(sph, gxyz, sumNS)
          import dbl, T_lateralGrid
          class(T_lateralGrid),   intent(in)    :: sph
          real(kind=dbl), target, intent(out)   :: gxyz(*)
          complex(kind=dbl),      intent(inout) :: sumNS(*)
        end subroutine grid_4_sub
        
        pure subroutine grid_8_sub(sph, gxyz, sumNS)
          import dbl, T_lateralGrid
          class(T_lateralGrid),   intent(in)    :: sph
          real(kind=dbl), target, intent(out)   :: gxyz(*)
          complex(kind=dbl),      intent(inout) :: sumNS(*)
        end subroutine grid_8_sub
        
        pure subroutine grid_16_sub(sph, gxyz, sumNS)
          import dbl, T_lateralGrid
          class(T_lateralGrid),   intent(in)    :: sph
          real(kind=dbl), target, intent(out)   :: gxyz(*)
          complex(kind=dbl),      intent(inout) :: sumNS(*)
        end subroutine grid_16_sub
      end interface
    end subroutine lege_transform_sub
    
    module pure subroutine vec2vec_jml_to_jml_sub(this, cjml, cab, ncab, cabpadding)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: ncab, cabpadding
      complex(kind=dbl),    intent(in)    :: cjml(*)
      complex(kind=dbl),    intent(inout) :: cab(ncab,*)
    end subroutine vec2vec_jml_to_jml_sub
    
    module pure subroutine scal2scal_jm_to_mj_sub(this, cjm, cab, ncab, cabpadding)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: ncab, cabpadding
      complex(kind=dbl),    intent(in)    :: cjm(*)
      complex(kind=dbl),    intent(inout) :: cab(ncab,*)
    end subroutine scal2scal_jm_to_mj_sub
    
    module pure subroutine vec2scal_jml_to_mj_sub(this, cab, ncab, cc)
      class(T_lateralGrid), intent(in)  :: this
      integer,              intent(in)  :: ncab
      complex(kind=dbl),    intent(in)  :: cab(ncab,*)
      complex(kind=dbl),    intent(out) :: cc(3,ncab,*)
    end subroutine vec2scal_jml_to_mj_sub
    
    module pure subroutine gradvec2vec_jmlk_to_jml_sub(this, ri, v, dv_r, cab, ncab, cabpadding)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: cabpadding, ncab
      real(kind=dbl),       intent(in)    :: ri
      complex(kind=dbl),    intent(in)    :: v(*), dv_r(*)
      complex(kind=dbl),    intent(inout) :: cab(ncab,*)
    end subroutine gradvec2vec_jmlk_to_jml_sub
    
    module pure subroutine devtens2scal_jml2_to_mj_sub(this, ctjml2, cr, ncr, crpadding)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: ncr, crpadding
      complex(kind=dbl),    intent(in)    :: ctjml2(*)
      complex(kind=dbl),    intent(inout) :: cr(ncr,*)
    end subroutine devtens2scal_jml2_to_mj_sub
    
    module pure subroutine scal2scal_mj_to_jm_sub(this, cr, ncr, crpadding, cjm, ncjm, cjmpadding)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: ncr, ncjm, crpadding, cjmpadding
      complex(kind=dbl),    intent(in)    :: cr(ncr,*)
      complex(kind=dbl),    intent(inout) :: cjm(ncjm,*)
    end subroutine scal2scal_mj_to_jm_sub
    
    module pure subroutine scal2vecscal_mj_to_jm_sub(this, cr, ncr, crpadding, cjm, ncjm, cjmpadding)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: ncr, crpadding, ncjm, cjmpadding
      complex(kind=dbl),    intent(inout) :: cr(ncr,*)
      complex(kind=dbl),    intent(inout) :: cjm(ncjm,*)
    end subroutine scal2vecscal_mj_to_jm_sub
    
    module pure subroutine scal2devtens_mj_to_jml2_sub(this, cr, ncr, crpadding, ctjml2)
      class(T_lateralGrid), intent(in) :: this
      integer,              intent(in)  :: ncr, crpadding
      complex(kind=dbl),    intent(in)  :: cr(ncr,*)
      complex(kind=dbl),    intent(out) :: ctjml2(*)
    end subroutine scal2devtens_mj_to_jml2_sub
    
    module subroutine vctol_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine vctol_sub
    
    module pure subroutine vcsum_sub(this, cajm, cbjm, cjm)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcsum_sub
    
    module pure subroutine vcst_sub(this, cajm, cbjml2, cjml2)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cajm(*), cbjml2(*)
      complex(kind=dbl),    intent(out) :: cjml2(*)
    end subroutine vcst_sub
    
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
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_2_vcsum_sub
    
    module pure subroutine grid_op_4_vcsum_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_4_vcsum_sub
    
    module pure subroutine grid_op_8_vcsum_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_8_vcsum_sub
    
    module pure subroutine grid_op_16_vcsum_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_16_vcsum_sub
    
    module pure subroutine grid_op_2_vcst_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_2_vcst_sub
    
    module pure subroutine grid_op_4_vcst_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_4_vcst_sub
    
    module pure subroutine grid_op_8_vcst_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_8_vcst_sub
    
    module pure subroutine grid_op_16_vcst_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_16_vcst_sub
    
    module pure subroutine grid_op_2_vcvv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_2_vcvv_sub
    
    module pure subroutine grid_op_4_vcvv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_4_vcvv_sub
    
    module pure subroutine grid_op_8_vcvv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_8_vcvv_sub
    
    module pure subroutine grid_op_16_vcvv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_16_vcvv_sub
    
    module pure subroutine grid_op_2_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_2_vcvgv_sub
    
    module pure subroutine grid_op_4_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_4_vcvgv_sub
    
    module pure subroutine grid_op_8_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_8_vcvgv_sub
    
    module pure subroutine grid_op_16_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_16_vcvgv_sub
    
    module pure subroutine grid_op_2_vcvv_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_2_vcvv_vcvgv_sub
    
    module pure subroutine grid_op_4_vcvv_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_4_vcvv_vcvgv_sub
    
    module pure subroutine grid_op_8_vcvv_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_8_vcvv_vcvgv_sub
    
    module pure subroutine grid_op_16_vcvv_vcvgv_sub(this, grid, sumNS)
      class(T_lateralGrid),   intent(in)    :: this
      real(kind=dbl), target, intent(out)   :: grid(*)
      complex(kind=dbl),      intent(inout) :: sumNS(*)
    end subroutine grid_op_16_vcvv_vcvgv_sub
  end interface
  
end module SphericalHarmonics
