module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                             :: jmax, nLege, nrma
    type(c_ptr)                         :: c_rw
    real(kind=dbl), pointer, contiguous :: rw(:,:)
    real(kind=dbl), allocatable         :: emj(:), fmj(:,:)
    
    contains
    
    procedure, public,  pass :: init_sub       => init_lege_sub
    procedure, private, pass :: roots_sub      => find_roots_sub
    procedure, private, pass :: coeffs_sub     => compute_coeffs_sub
    procedure, public,  pass :: deallocate_sub => deallocate_lege_sub
    
    procedure, public,  pass :: alloc_rscal_sub => allocate_rscalars_sub
    procedure, public,  pass :: index_bwd_sub   => c2r_mj_to_mj_sub
    procedure, public,  pass :: index_fwd_sub   => r2c_mj_to_mj_sub
    
    procedure, public,  pass :: bwd_legesum_sub
    procedure, public,  pass :: fwd_legesum_sub
    
  end type T_legep
  
  interface
    module subroutine init_lege_sub(this, jmax, nLege, wfac)
      class(T_legep), intent(inout) :: this
      integer,        intent(in)    :: jmax, nLege
      real(kind=dbl), intent(in)    :: wfac
    end subroutine init_lege_sub
    
    module  subroutine deallocate_lege_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine deallocate_lege_sub
    
    module subroutine find_roots_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine find_roots_sub
    
    module  subroutine compute_coeffs_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine compute_coeffs_sub
    
    module  subroutine allocate_rscalars_sub(this, ns, rscal)
      class(T_legep),              intent(in)  :: this
      integer,                     intent(in)  :: ns
      real(kind=dbl), allocatable, intent(out) :: rscal(:)
    end subroutine allocate_rscalars_sub
    
    module  subroutine c2r_mj_to_mj_sub(this, ncab, cab, rcab)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: ncab
      complex(kind=dbl), intent(in)  :: cab(ncab,*)
      real(kind=dbl),    intent(out) :: rcab(2,ncab,2,this%nrma)
    end subroutine c2r_mj_to_mj_sub
    
    module  subroutine r2c_mj_to_mj_sub(this, ncab, cab, rcab)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: ncab
      real(kind=dbl),    intent(in)  :: rcab(2,ncab,2,this%nrma)
      complex(kind=dbl), intent(out) :: cab(ncab,*)
    end subroutine r2c_mj_to_mj_sub
    
    module  subroutine bwd_legesum_sub(this, nb, cc, sumN, sumS, cosx, sinx, cosx2, pmm, pmj2, pmj1, pmj, swork)
      class(T_legep), intent(in)  :: this
      integer,        intent(in)  :: nb
      real(kind=dbl), intent(in)  :: cosx(step), sinx(step), cosx2(step)
      real(kind=dbl), intent(out) :: pmm(step), pmj2(step), pmj1(step), pmj(step), swork(4*nb*step)
      real(kind=dbl), intent(out) :: sumN(2*nb*step,0:this%jmax), sumS(2*nb*step,0:this%jmax)
      real(kind=dbl), intent(in)  :: cc(4*nb,this%nrma)
    end subroutine bwd_legesum_sub
    
    module  subroutine fwd_legesum_sub(this, nf, sumN, sumS, cr, cosx, sinx, cosx2, weight, pmm, pmj2, pmj1, pmj, swork)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: nf
      real(kind=dbl), intent(in)    :: cosx(step), sinx(step), cosx2(step), weight(step)
      real(kind=dbl), intent(out)   :: pmm(step), pmj2(step), pmj1(step), pmj(step), swork(step)
      real(kind=dbl), intent(in)    :: sumN(2*nf*step,0:this%jmax), sumS(2*nf*step,0:this%jmax)
      real(kind=dbl), intent(inout) :: cr(4*nf,this%nrma)
    end subroutine fwd_legesum_sub
  end interface
  
  interface
    module pure subroutine mmset_sub(ma, cff, cosx, sinx, pmm, pmj2, pmj1, pmj)
      integer,        intent(in)    :: ma
      real(kind=dbl), intent(in)    :: cff, cosx(step), sinx(step)
      real(kind=dbl), intent(inout) :: pmm(step)
      real(kind=dbl), intent(out)   :: pmj2(step), pmj1(step), pmj(step)
    end subroutine mmset_sub
    
    module pure subroutine mjrec_sub(cff, cosx2, pmj2, pmj1, pmj)
      real(kind=dbl), intent(in)    :: cff(3), cosx2(step)
      real(kind=dbl), intent(inout) :: pmj2(step), pmj1(step), pmj(step)
    end subroutine mjrec_sub
    
    module pure subroutine bwd_sum_sub(n, pmj, cc, swork)
      integer,        intent(in)  :: n
      real(kind=dbl), intent(in)  :: pmj(step)
      real(kind=dbl), intent(in)  :: cc(n)
      real(kind=dbl), intent(out) :: swork(step,n)
    end subroutine bwd_sum_sub
    
    module pure subroutine bwd_shuffle_sub(n, cosx, swork, sumN, sumS)
      integer,        intent(in)    :: n
      real(kind=dbl), intent(in)    :: cosx(step)
      real(kind=dbl), intent(inout) :: swork(step,2,n,2)
      real(kind=dbl), intent(out)   :: sumN(step,n,2), sumS(step,n,2)
    end subroutine bwd_shuffle_sub
    
    module pure subroutine fwd_sum_sub(n, pmj, swork, cr)
      integer,        intent(in)    :: n
      real(kind=dbl), intent(in)    :: pmj(step)
      real(kind=dbl), intent(in)    :: swork(step,n)
      real(kind=dbl), intent(inout) :: cr(n)
    end subroutine fwd_sum_sub
    
    module pure subroutine fwd_shuffle_sub(n, w, cosx, sumN, sumS, swork)
      integer,        intent(in)  :: n
      real(kind=dbl), intent(in)  :: w(step), cosx(step), sumN(step,n,2), sumS(step,n,2)
      real(kind=dbl), intent(out) :: swork(step,2,n,2)
    end subroutine fwd_shuffle_sub
  end interface
  
end module lege_poly