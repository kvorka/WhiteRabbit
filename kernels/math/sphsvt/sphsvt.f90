module sphsvt
  use clebsch_gordan
  implicit none
  
  type, public :: T_sphsvt
    integer :: jmax, jmax1, jmax2, jmax3, jms, jms1, jms2, jmv, jmv1, jmv2
    
    contains
    
    procedure :: init_sub => init_sphsvt_sub
    procedure :: scal2scal_jm_to_mj_sub, scal2scal_mj_to_jm_sub, vec2vec_jml_to_jml_sub, vec2scal_jml_to_mj_sub, &
               & scal2vecscal_mj_to_jm_sub, gradvec2vec_jmlk_to_jml_sub, devtens2scal_jml2_to_mj_sub, &
               & scal2devtens_mj_to_jml2_sub, scal2vec_mj_to_jml_sub
    
  end type T_sphsvt
  
  interface
    module pure subroutine init_sphsvt_sub(this, jmax)
      class(T_sphsvt), intent(inout) :: this
      integer,         intent(in)    :: jmax
    end subroutine init_sphsvt_sub
    
    module pure subroutine vec2vec_jml_to_jml_sub(this, cjml, cab, ncab, cabpadding)
      class(T_sphsvt),   intent(in)    :: this
      integer,           intent(in)    :: ncab, cabpadding
      complex(kind=dbl), intent(in)    :: cjml(*)
      complex(kind=dbl), intent(inout) :: cab(ncab,*)
    end subroutine vec2vec_jml_to_jml_sub
    
    module pure subroutine scal2scal_jm_to_mj_sub(this, cjm, cab, ncab, cabpadding)
      class(T_sphsvt),   intent(in)    :: this
      integer,           intent(in)    :: ncab, cabpadding
      complex(kind=dbl), intent(in)    :: cjm(*)
      complex(kind=dbl), intent(inout) :: cab(ncab,*)
    end subroutine scal2scal_jm_to_mj_sub
    
    module pure subroutine vec2scal_jml_to_mj_sub(this, cab, ncab, cc, ncc, ccpadding)
      class(T_sphsvt),   intent(in)    :: this
      integer,           intent(in)    :: ncab, ncc, ccpadding
      complex(kind=dbl), intent(in)    :: cab(ncab,*)
      complex(kind=dbl), intent(inout) :: cc(ncc,*)
    end subroutine vec2scal_jml_to_mj_sub
    
    module pure subroutine gradvec2vec_jmlk_to_jml_sub(this, ri, v, dv_r, cab, ncab, cabpadding)
      class(T_sphsvt),   intent(in)    :: this
      integer,           intent(in)    :: cabpadding, ncab
      real(kind=dbl),    intent(in)    :: ri
      complex(kind=dbl), intent(in)    :: v(*), dv_r(*)
      complex(kind=dbl), intent(inout) :: cab(ncab,*)
    end subroutine gradvec2vec_jmlk_to_jml_sub
    
    module pure subroutine devtens2scal_jml2_to_mj_sub(this, ctjml2, cr, ncr, crpadding)
      class(T_sphsvt),   intent(in)    :: this
      integer,           intent(in)    :: ncr, crpadding
      complex(kind=dbl), intent(in)    :: ctjml2(*)
      complex(kind=dbl), intent(inout) :: cr(ncr,*)
    end subroutine devtens2scal_jml2_to_mj_sub
    
    module pure subroutine scal2scal_mj_to_jm_sub(this, cr, ncr, crpadding, cjm, ncjm, cjmpadding)
      class(T_sphsvt),   intent(in)    :: this
      integer,           intent(in)    :: ncr, ncjm, crpadding, cjmpadding
      complex(kind=dbl), intent(in)    :: cr(ncr,*)
      complex(kind=dbl), intent(inout) :: cjm(ncjm,*)
    end subroutine scal2scal_mj_to_jm_sub
    
    module pure subroutine scal2vecscal_mj_to_jm_sub(this, cr, ncr, crpadding, cjm, ncjm, cjmpadding)
      class(T_sphsvt),   intent(in)    :: this
      integer,           intent(in)    :: ncr, crpadding, ncjm, cjmpadding
      complex(kind=dbl), intent(inout) :: cr(ncr,*)
      complex(kind=dbl), intent(inout) :: cjm(ncjm,*)
    end subroutine scal2vecscal_mj_to_jm_sub
    
    module pure subroutine scal2vec_mj_to_jml_sub(this, cr, ncr, crpadding, cjml, ncjml, cjmlpadding)
      class(T_sphsvt),   intent(in)    :: this
      integer,           intent(in)    :: ncr, crpadding, ncjml, cjmlpadding
      complex(kind=dbl), intent(inout) :: cr(ncr,*)
      complex(kind=dbl), intent(inout) :: cjml(ncjml,*)
    end subroutine scal2vec_mj_to_jml_sub
    
    module pure subroutine scal2devtens_mj_to_jml2_sub(this, cr, ncr, crpadding, ctjml2)
      class(T_sphsvt),   intent(in) :: this
      integer,           intent(in)  :: ncr, crpadding
      complex(kind=dbl), intent(in)  :: cr(ncr,*)
      complex(kind=dbl), intent(out) :: ctjml2(*)
    end subroutine scal2devtens_mj_to_jml2_sub
  end interface
  
end module sphsvt