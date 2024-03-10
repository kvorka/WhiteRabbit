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
    integer,                     private :: jmax, jmax1, jmax2, jmax3, jms, jms1, jms2, jmv, jmv1, jmt, nLegendre, nFourier
    integer,        allocatable, private :: maxm(:)
    real(kind=dbl),              private :: tolm, scale
    real(kind=dbl), allocatable, private :: roots(:), fftLege(:), amjrr(:), bmjrr(:), pmm(:,:)
    type(T_fft),                 private :: fourtrans
    
    contains
    
    !Allocation+initialization and cleaning
    procedure :: init_sub       => init_harmonics_sub
    procedure :: deallocate_sub => deallocate_harmonics_sub
    
    !Vector transforms, scalar reindexing
    procedure :: scal2scal_jm_to_mj_sub, scal2scal_mj_to_jm_sub, vec2vec_jml_to_jml_sub, vec2scal_jml_to_mj_sub, &
               & scal2vecscal_mj_to_jm_sub, gradvec2vec_jmlk_to_jml_sub, devtens2scal_jml2_to_mj_sub, scal2devtens_mj_to_jml2_sub
    
    !Generic scalar transforms and needed routines
    procedure :: get_maxm_fn, rescale_sub, lege_transform_sub
    
    !Precision setting and transforms
    procedure :: vctol_sub, vcsum_sub, vcst_sub, vcvv_sub, vcvgv_sub, vcvv_vcvgv_sub
    
  end type T_lateralGrid
  
  interface
    module subroutine init_harmonics_sub(this, jmax)
      class(T_lateralGrid), intent(inout) :: this
      integer,              intent(in)    :: jmax
    end subroutine init_harmonics_sub
    
    module pure subroutine deallocate_harmonics_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_harmonics_sub
    
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
    
    module pure subroutine rescale_sub( this, arr, length )
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: length
      complex(kind=dbl),    intent(inout) :: arr(*)
    end subroutine rescale_sub
    
    module pure integer function get_maxm_fn(this, i, i2)
      class(T_lateralGrid), intent(in) :: this
      integer,              intent(in) :: i, i2
    end function get_maxm_fn
    
    module pure subroutine lege_transform_sub(this, nforw, nback, cc, cr, grid_sub)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: nforw, nback
      complex(kind=dbl),    intent(in)    :: cc(nback,*)
      complex(kind=dbl),    intent(inout) :: cr(nforw,*)
      
      interface
        module pure subroutine grid_sub(nfour, nstep, gxyz)
          integer,                intent(in)    :: nfour, nstep
          real(kind=dbl), target, intent(inout) :: gxyz(*)
        end subroutine grid_sub
      end interface
    end subroutine lege_transform_sub
    
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
  end interface
  
end module SphericalHarmonics
