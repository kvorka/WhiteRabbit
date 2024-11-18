module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                     :: jmax, nLege
    real(kind=dbl), allocatable :: abmj(:,:), rw(:,:)
    
    contains
    
    procedure, public,  pass :: init_sub       => init_lege_sub
    procedure, private, pass :: roots_sub      => find_roots_sub
    procedure, private, pass :: coeffs_sub     => compute_coeffs_sub
    procedure, public,  pass :: deallocate_sub => deallocate_lege_sub
    
    procedure, private, pass :: mmset_4_sub, recursion_4_sub, recursion2_4_sub
    procedure, private, pass :: backward_sum_4_sub, backward_sum2_4_sub, backward_rcb_4_sub
    procedure, private, pass :: forward_sum_4_sub, forward_sum2_4_sub, forward_rcb_4_sub
    procedure, public,  pass :: backward_legesum_4_sub, forward_legesum_4_sub
    
    procedure, private, pass :: mmset_8_sub, recursion_8_sub, recursion2_8_sub
    procedure, private, pass :: backward_sum_8_sub, backward_sum2_8_sub, backward_rcb_8_sub
    procedure, private, pass :: forward_sum_8_sub, forward_sum2_8_sub, forward_rcb_8_sub
    procedure, public,  pass :: backward_legesum_8_sub, forward_legesum_8_sub
    
    procedure, private, pass :: mmset_16_sub, recursion_16_sub, recursion2_16_sub
    procedure, private, pass :: backward_sum_16_sub, backward_sum2_16_sub, backward_rcb_16_sub
    procedure, private, pass :: forward_sum_16_sub, forward_sum2_16_sub, forward_rcb_16_sub
    procedure, public,  pass :: backward_legesum_16_sub, forward_legesum_16_sub
    
  end type T_legep
  
  interface
    module subroutine init_lege_sub(this, jmax, nLege, wfac)
      class(T_legep), intent(inout) :: this
      integer,        intent(in)    :: jmax, nLege
      real(kind=dbl), intent(in)    :: wfac
    end subroutine init_lege_sub
    
    module pure subroutine deallocate_lege_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine deallocate_lege_sub
    
    module subroutine find_roots_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine find_roots_sub
    
    module pure subroutine compute_coeffs_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine compute_coeffs_sub
    
    module pure subroutine mmset_4_sub(this, mj, sinx, pmm, pmj)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: sinx(4)
      real(kind=dbl), intent(inout) :: pmm(4)
      real(kind=dbl), intent(out)   :: pmj(4,3)
    end subroutine mmset_4_sub
    
    module pure subroutine recursion_4_sub(this, mj, cosx, pmj)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: cosx(4)
      real(kind=dbl), intent(inout) :: pmj(4,3)
    end subroutine recursion_4_sub
    
    module pure subroutine recursion2_4_sub(this, mj, cosx, pmj)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: cosx(4)
      real(kind=dbl), intent(inout) :: pmj(4,3)
    end subroutine recursion2_4_sub
    
    module pure subroutine backward_sum_4_sub(this, nb, legep, cc, legesum)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nb
      real(kind=dbl),    intent(in)    :: legep(4)
      complex(kind=dbl), intent(in)    :: cc(nb)
      complex(kind=dbl), intent(inout) :: legesum(4,nb)
    end subroutine backward_sum_4_sub
    
    module pure subroutine backward_sum2_4_sub(this, nb, legep, cc, legesum)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nb
      real(kind=dbl),    intent(in)    :: legep(4,2)
      complex(kind=dbl), intent(in)    :: cc(nb,2)
      complex(kind=dbl), intent(inout) :: legesum(4,nb,2)
    end subroutine backward_sum2_4_sub
    
    module pure subroutine backward_rcb_4_sub(this, nb, swork, sumN, sumS)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nb
      complex(kind=dbl), intent(in)  :: swork(4,nb,2)
      real(kind=dbl),    intent(out) :: sumN(4,nb,2), sumS(4,nb,2)
    end subroutine backward_rcb_4_sub
    
    module pure subroutine forward_sum_4_sub(this, nf, legep, legesum, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: legep(4)
      complex(kind=dbl), intent(in)    :: legesum(4,nf)
      complex(kind=dbl), intent(inout) :: cr(nf)
    end subroutine forward_sum_4_sub
    
    module pure subroutine forward_sum2_4_sub(this, nf, legep, legesum, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: legep(4,2)
      complex(kind=dbl), intent(in)    :: legesum(4,nf,2)
      complex(kind=dbl), intent(inout) :: cr(nf,2)
    end subroutine forward_sum2_4_sub
    
    module pure subroutine forward_rcb_4_sub(this, nf, w, sumN, sumS, swork)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nf
      real(kind=dbl),    intent(in)  :: w(4)
      real(kind=dbl),    intent(in)  :: sumN(4,nf,2), sumS(4,nf,2)
      complex(kind=dbl), intent(out) :: swork(4,nf,2)
    end subroutine forward_rcb_4_sub
    
    module pure subroutine backward_legesum_4_sub(this, nb, cc, sumN, sumS, rw)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nb
      complex(kind=dbl), intent(in)  :: cc(nb,*)
      real(kind=dbl),    intent(in)  :: rw(4,2)
      real(kind=dbl),    intent(out) :: sumN(8*nb,0:this%jmax), sumS(8*nb,0:this%jmax)
    end subroutine backward_legesum_4_sub
    
    module pure subroutine forward_legesum_4_sub(this, nf, sumN, sumS, cr, rw)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: sumN(8*nf,0:this%jmax), sumS(8*nf,0:this%jmax), rw(4,3)
      complex(kind=dbl), intent(inout) :: cr(nf,*)
    end subroutine forward_legesum_4_sub
    
    module pure subroutine mmset_8_sub(this, mj, sinx, pmm, pmj)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: sinx(8)
      real(kind=dbl), intent(inout) :: pmm(8)
      real(kind=dbl), intent(out)   :: pmj(8,3)
    end subroutine mmset_8_sub
    
    module pure subroutine recursion_8_sub(this, mj, cosx, pmj)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: cosx(8)
      real(kind=dbl), intent(inout) :: pmj(8,3)
    end subroutine recursion_8_sub
    
    module pure subroutine recursion2_8_sub(this, mj, cosx, pmj)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: cosx(8)
      real(kind=dbl), intent(inout) :: pmj(8,3)
    end subroutine recursion2_8_sub
    
    module pure subroutine backward_sum_8_sub(this, nb, legep, cc, legesum)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nb
      real(kind=dbl),    intent(in)    :: legep(8)
      complex(kind=dbl), intent(in)    :: cc(nb)
      complex(kind=dbl), intent(inout) :: legesum(8,nb)
    end subroutine backward_sum_8_sub
    
    module pure subroutine backward_sum2_8_sub(this, nb, legep, cc, legesum)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nb
      real(kind=dbl),    intent(in)    :: legep(8,2)
      complex(kind=dbl), intent(in)    :: cc(nb,2)
      complex(kind=dbl), intent(inout) :: legesum(8,nb,2)
    end subroutine backward_sum2_8_sub
    
    module pure subroutine backward_rcb_8_sub(this, nb, swork, sumN, sumS)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nb
      complex(kind=dbl), intent(in)  :: swork(8,nb,2)
      real(kind=dbl),    intent(out) :: sumN(8,nb,2), sumS(8,nb,2)
    end subroutine backward_rcb_8_sub
    
    module pure subroutine forward_sum_8_sub(this, nf, legep, legesum, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: legep(8)
      complex(kind=dbl), intent(in)    :: legesum(8,nf)
      complex(kind=dbl), intent(inout) :: cr(nf)
    end subroutine forward_sum_8_sub
    
    module pure subroutine forward_sum2_8_sub(this, nf, legep, legesum, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: legep(8,2)
      complex(kind=dbl), intent(in)    :: legesum(8,nf,2)
      complex(kind=dbl), intent(inout) :: cr(nf,2)
    end subroutine forward_sum2_8_sub
    
    module pure subroutine forward_rcb_8_sub(this, nf, w, sumN, sumS, swork)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nf
      real(kind=dbl),    intent(in)  :: w(8)
      real(kind=dbl),    intent(in)  :: sumN(8,nf,2), sumS(8,nf,2)
      complex(kind=dbl), intent(out) :: swork(8,nf,2)
    end subroutine forward_rcb_8_sub
    
    module pure subroutine backward_legesum_8_sub(this, nb, cc, sumN, sumS, rw)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nb
      real(kind=dbl),    intent(in)  :: rw(8,2)
      real(kind=dbl),    intent(out) :: sumN(16*nb,0:this%jmax), sumS(16*nb,0:this%jmax)
      complex(kind=dbl), intent(in)  :: cc(nb,*)
    end subroutine backward_legesum_8_sub
    
    module pure subroutine forward_legesum_8_sub(this, nf, sumN, sumS, cr, rw)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: sumN(16*nf,0:this%jmax), sumS(16*nf,0:this%jmax), rw(8,3)
      complex(kind=dbl), intent(inout) :: cr(nf,*)
    end subroutine forward_legesum_8_sub
    
    module pure subroutine mmset_16_sub(this, mj, sinx, pmm, pmj)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: sinx(16)
      real(kind=dbl), intent(inout) :: pmm(16)
      real(kind=dbl), intent(out)   :: pmj(16,3)
    end subroutine mmset_16_sub
    
    module pure subroutine recursion_16_sub(this, mj, cosx, pmj)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: cosx(16)
      real(kind=dbl), intent(inout) :: pmj(16,3)
    end subroutine recursion_16_sub
    
    module pure subroutine recursion2_16_sub(this, mj, cosx, pmj)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: cosx(16)
      real(kind=dbl), intent(inout) :: pmj(16,3)
    end subroutine recursion2_16_sub
    
    module pure subroutine backward_sum_16_sub(this, nb, legep, cc, legesum)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nb
      real(kind=dbl),    intent(in)    :: legep(16)
      complex(kind=dbl), intent(in)    :: cc(nb)
      complex(kind=dbl), intent(inout) :: legesum(16,nb)
    end subroutine backward_sum_16_sub
    
    module pure subroutine backward_sum2_16_sub(this, nb, legep, cc, legesum)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nb
      real(kind=dbl),    intent(in)    :: legep(16,2)
      complex(kind=dbl), intent(in)    :: cc(nb,2)
      complex(kind=dbl), intent(inout) :: legesum(16,nb,2)
    end subroutine backward_sum2_16_sub
    
    module pure subroutine backward_rcb_16_sub(this, nb, swork, sumN, sumS)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nb
      complex(kind=dbl), intent(in)  :: swork(16,nb,2)
      real(kind=dbl),    intent(out) :: sumN(16,nb,2), sumS(16,nb,2)
    end subroutine backward_rcb_16_sub
    
    module pure subroutine forward_sum_16_sub(this, nf, legep, legesum, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: legep(16)
      complex(kind=dbl), intent(in)    :: legesum(16,nf)
      complex(kind=dbl), intent(inout) :: cr(nf)
    end subroutine forward_sum_16_sub
    
    module pure subroutine forward_sum2_16_sub(this, nf, legep, legesum, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: legep(16,2)
      complex(kind=dbl), intent(in)    :: legesum(16,nf,2)
      complex(kind=dbl), intent(inout) :: cr(nf,2)
    end subroutine forward_sum2_16_sub
    
    module pure subroutine forward_rcb_16_sub(this, nf, w, sumN, sumS, swork)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nf
      real(kind=dbl),    intent(in)  :: w(16)
      real(kind=dbl),    intent(in)  :: sumN(16,nf,2), sumS(16,nf,2)
      complex(kind=dbl), intent(out) :: swork(16,nf,2)
    end subroutine forward_rcb_16_sub
    
    module pure subroutine backward_legesum_16_sub(this, nb, cc, sumN, sumS, rw)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nb
      real(kind=dbl),    intent(in)  :: rw(16,2)
      real(kind=dbl),    intent(out) :: sumN(32*nb,0:this%jmax), sumS(32*nb,0:this%jmax)
      complex(kind=dbl), intent(in)  :: cc(nb,*)
    end subroutine backward_legesum_16_sub
    
    module pure subroutine forward_legesum_16_sub(this, nf, sumN, sumS, cr, rw)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: sumN(32*nf,0:this%jmax), sumS(32*nf,0:this%jmax), rw(16,3)
      complex(kind=dbl), intent(inout) :: cr(nf,*)
    end subroutine forward_legesum_16_sub
  end interface
  
end module lege_poly