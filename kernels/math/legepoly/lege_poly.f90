module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                     :: jmax, nLege, nrma
    real(kind=dbl), pointer     :: rw(:,:)
    real(kind=dbl), allocatable :: emj(:), fmj(:,:)
    type(c_ptr)                 :: c_rw
    
    contains
    
    procedure, public,  pass :: init_sub       => init_lege_sub
    procedure, private, pass :: roots_sub      => find_roots_sub
    procedure, private, pass :: coeffs_sub     => compute_coeffs_sub
    procedure, public,  pass :: deallocate_sub => deallocate_lege_sub
    
    procedure, public,  pass :: alloc_cscal_sub => allocate_cscalars_sub
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
    
    module  subroutine allocate_cscalars_sub(this, ns, cscal)
      class(T_legep),                 intent(in)  :: this
      integer,                        intent(in)  :: ns
      complex(kind=dbl), allocatable, intent(out) :: cscal(:)
    end subroutine allocate_cscalars_sub
    
    module  subroutine c2r_mj_to_mj_sub(this, ncab, cab, rcab)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: ncab
      complex(kind=dbl), intent(in)  :: cab(ncab,*)
      complex(kind=dbl), intent(out) :: rcab(ncab,2,this%nrma)
    end subroutine c2r_mj_to_mj_sub
    
    module  subroutine r2c_mj_to_mj_sub(this, ncab, cab, rcab)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: ncab
      complex(kind=dbl), intent(in)  :: rcab(ncab,2,this%nrma)
      complex(kind=dbl), intent(out) :: cab(ncab,*)
    end subroutine r2c_mj_to_mj_sub
    
    module  subroutine bwd_legesum_sub(this, nb, cc, sumN, sumS, cosx, sinx, cosx2, pmm, pmj2, pmj1, pmj, swork)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nb
      real(kind=dbl),    intent(in)  :: cosx(*), sinx(*), cosx2(*)
      real(kind=dbl),    intent(out) :: pmm(*), pmj2(*), pmj1(*), pmj(*), swork(*)
      real(kind=dbl),    intent(out) :: sumN(2*nb*step,0:*), sumS(2*nb*step,0:*)
      complex(kind=dbl), intent(in)  :: cc(2*nb,*)
    end subroutine bwd_legesum_sub
    
    module  subroutine fwd_legesum_sub(this, nf, sumN, sumS, cr, cosx, sinx, cosx2, weight, pmm, pmj2, pmj1, pmj, swork)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: cosx(step), sinx(step), cosx2(step), weight(step)
      real(kind=dbl),    intent(out)   :: pmm(step), pmj2(step), pmj1(step), pmj(step), swork(4*nf*step)
      real(kind=dbl),    intent(in)    :: sumN(2*nf*step,0:this%jmax), sumS(2*nf*step,0:this%jmax)
      complex(kind=dbl), intent(inout) :: cr(2*nf,this%nrma)
    end subroutine fwd_legesum_sub
  end interface
  
  interface
    subroutine mmset_sub(n1, n2, cff, cosx, sinx, pmm, pmj2, pmj1, pmj) bind(C, name='mmset')
      import                        :: dbl
      integer,        intent(in)    :: n1, n2
      real(kind=dbl), intent(in)    :: cff, cosx(*), sinx(*)
      real(kind=dbl), intent(inout) :: pmm(*)
      real(kind=dbl), intent(out)   :: pmj2(*), pmj1(*), pmj(*)
    end subroutine mmset_sub
    
    subroutine mjrec_sub(n1, cff, cosx2, pmj2, pmj1, pmj) bind(C, name='mjrec')
      import                        :: dbl
      integer,        intent(in)    :: n1
      real(kind=dbl), intent(in)    :: cff(*), cosx2(*)
      real(kind=dbl), intent(inout) :: pmj2(*), pmj1(*), pmj(*)
    end subroutine mjrec_sub
    
    subroutine bwd_sum_sub(n1, n2, pmj, cc, swork) bind(C, name='bwd_sum')
      import                           :: dbl
      integer,           intent(in)    :: n1, n2
      real(kind=dbl),    intent(in)    :: pmj(*)
      complex(kind=dbl), intent(in)    :: cc(*)
      real(kind=dbl),    intent(inout) :: swork(*)
    end subroutine bwd_sum_sub
    
    subroutine bwd_shuffle_sub(n1, n2, cosx, swork, sumN, sumS) bind(C, name='bwd_shuffle')
      import                        :: dbl
      integer,        intent(in)    :: n1, n2
      real(kind=dbl), intent(in)    :: cosx(*)
      real(kind=dbl), intent(inout) :: swork(*)
      real(kind=dbl), intent(out)   :: sumN(*), sumS(*)
    end subroutine bwd_shuffle_sub
    
    subroutine fwd_sum_sub(n1, n2, pmj, swork, cr) bind(C, name='fwd_sum')
      import                           :: dbl
      integer,           intent(in)    :: n1, n2
      real(kind=dbl),    intent(in)    :: pmj(*), swork(*)
      complex(kind=dbl), intent(inout) :: cr(*)
    end subroutine fwd_sum_sub
    
    subroutine fwd_shuffle_sub(n1, n2, wght, cosx, sumN, sumS, swork) bind(C, name='fwd_shuffle')
      import                      :: dbl
      integer,        intent(in)  :: n1, n2
      real(kind=dbl), intent(in)  :: wght(*), cosx(*), sumN(*), sumS(*)
      real(kind=dbl), intent(out) :: swork(*)
    end subroutine fwd_shuffle_sub
  end interface
  
end module lege_poly