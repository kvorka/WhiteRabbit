module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                     :: jmax, nLege
    real(kind=dbl), allocatable :: amj(:), bmj(:), cmm(:), rootsweights(:,:)
    
    contains
    
    procedure, public,  pass :: init_sub       => init_lege_sub
    procedure, private, pass :: roots_sub      => find_roots_sub
    procedure, private, pass :: coeffs_sub     => compute_coeffs_sub
    procedure, public,  pass :: deallocate_sub => deallocate_lege_sub
    
    procedure, private, pass :: mmset_4_sub, recursion_4_sub
    procedure, private, pass :: backward_sum_4_sub, backward_rcb_4_sub
    procedure, private, pass :: forward_sum_4_sub, forward_rcb_4_sub
    procedure, public,  pass :: backward_legesum_4_sub, forward_legesum_4_sub
    
    procedure, private, pass :: mmset_8_sub, recursion_8_sub
    procedure, private, pass :: backward_sum_8_sub, backward_rcb_8_sub
    procedure, private, pass :: forward_sum_8_sub, forward_rcb_8_sub
    procedure, public,  pass :: backward_legesum_8_sub, forward_legesum_8_sub
    
    procedure, private, pass :: mmset_16_sub, recursion_16_sub
    procedure, private, pass :: backward_sum_16_sub, backward_rcb_16_sub
    procedure, private, pass :: forward_sum_16_sub, forward_rcb_16_sub
    procedure, public,  pass :: backward_legesum_16_sub, forward_legesum_16_sub
    
  end type T_legep
  
  interface
    module subroutine init_lege_sub(this, jmax, nLege)
      class(T_legep), intent(inout) :: this
      integer,        intent(in)    :: jmax, nLege
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
    
    module pure subroutine mmset_4_sub(this, m, sinx, pmm, pmj2, pmj1, pmj0)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: m
      real(kind=dbl), intent(in)    :: sinx(4)
      real(kind=dbl), intent(inout) :: pmm(4)
      real(kind=dbl), intent(out)   :: pmj2(4), pmj1(4), pmj0(4)
    end subroutine mmset_4_sub
    
    module pure subroutine recursion_4_sub(this, mj, cosx, pmj2, pmj1, pmj0)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: cosx(4)
      real(kind=dbl), intent(inout) :: pmj2(4), pmj1(4), pmj0(4)
    end subroutine recursion_4_sub
    
    module pure subroutine backward_sum_4_sub(this, nb, legep, cc, legesum)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nb
      real(kind=dbl),    intent(in)    :: legep(4)
      complex(kind=dbl), intent(in)    :: cc(nb)
      complex(kind=dbl), intent(inout) :: legesum(4,nb)
    end subroutine backward_sum_4_sub
    
    module pure subroutine backward_rcb_4_sub(this, nb, sumsym, sumasym, sumN, sumS)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nb
      complex(kind=dbl), intent(in)  :: sumsym(4,nb), sumasym(4,nb)
      real(kind=dbl),    intent(out) :: sumN(nb,4,2), sumS(nb,4,2)
    end subroutine backward_rcb_4_sub
    
    module pure subroutine forward_sum_4_sub(this, nf, legep, legesum, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: legep(4)
      complex(kind=dbl), intent(in)    :: legesum(4,nf)
      complex(kind=dbl), intent(inout) :: cr(nf)
    end subroutine forward_sum_4_sub
    
    module pure subroutine forward_rcb_4_sub(this, nf, w, sumN, sumS, sumsym, sumasym)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nf
      real(kind=dbl),    intent(in)  :: w(4)
      real(kind=dbl),    intent(in)  :: sumN(nf,4,2), sumS(nf,4,2)
      complex(kind=dbl), intent(out) :: sumsym(4,nf), sumasym(4,nf)
    end subroutine forward_rcb_4_sub
    
    module pure subroutine backward_legesum_4_sub(this, it, nb, cc, sumN, sumS)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: it, nb
      complex(kind=dbl), intent(in)  :: cc(nb,*)
      real(kind=dbl),    intent(out) :: sumN(nb,4,2,0:this%jmax), sumS(nb,4,2,0:this%jmax)
    end subroutine backward_legesum_4_sub
    
    module pure subroutine forward_legesum_4_sub(this, it, nf, sumN, sumS, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: it, nf
      real(kind=dbl),    intent(in)    :: sumN(nf,4,2,0:this%jmax), sumS(nf,4,2,0:this%jmax)
      complex(kind=dbl), intent(inout) :: cr(nf,*)
    end subroutine forward_legesum_4_sub
    
    module pure subroutine mmset_8_sub(this, m, sinx, pmm, pmj2, pmj1, pmj0)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: m
      real(kind=dbl), intent(in)    :: sinx(8)
      real(kind=dbl), intent(inout) :: pmm(8)
      real(kind=dbl), intent(out)   :: pmj2(8), pmj1(8), pmj0(8)
    end subroutine mmset_8_sub
    
    module pure subroutine recursion_8_sub(this, mj, cosx, pmj2, pmj1, pmj0)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: cosx(8)
      real(kind=dbl), intent(inout) :: pmj2(8), pmj1(8), pmj0(8)
    end subroutine recursion_8_sub
    
    module pure subroutine backward_sum_8_sub(this, nb, legep, cc, legesum)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nb
      real(kind=dbl),    intent(in)    :: legep(8)
      complex(kind=dbl), intent(in)    :: cc(nb)
      complex(kind=dbl), intent(inout) :: legesum(8,nb)
    end subroutine backward_sum_8_sub
    
    module pure subroutine backward_rcb_8_sub(this, nb, sumsym, sumasym, sumN, sumS)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nb
      complex(kind=dbl), intent(in)  :: sumsym(8,nb), sumasym(8,nb)
      real(kind=dbl),    intent(out) :: sumN(nb,8,2), sumS(nb,8,2)
    end subroutine backward_rcb_8_sub
    
    module pure subroutine forward_sum_8_sub(this, nf, legep, legesum, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: legep(8)
      complex(kind=dbl), intent(in)    :: legesum(8,nf)
      complex(kind=dbl), intent(inout) :: cr(nf)
    end subroutine forward_sum_8_sub
    
    module pure subroutine forward_rcb_8_sub(this, nf, w, sumN, sumS, sumsym, sumasym)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nf
      real(kind=dbl),    intent(in)  :: w(8)
      real(kind=dbl),    intent(in)  :: sumN(nf,8,2), sumS(nf,8,2)
      complex(kind=dbl), intent(out) :: sumsym(8,nf), sumasym(8,nf)
    end subroutine forward_rcb_8_sub
    
    module pure subroutine backward_legesum_8_sub(this, it, nb, cc, sumN, sumS)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: it, nb
      complex(kind=dbl), intent(in)  :: cc(nb,*)
      real(kind=dbl),    intent(out) :: sumN(nb,8,2,0:this%jmax), sumS(nb,8,2,0:this%jmax)
    end subroutine backward_legesum_8_sub
    
    module pure subroutine forward_legesum_8_sub(this, it, nf, sumN, sumS, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: it, nf
      real(kind=dbl),    intent(in)    :: sumN(nf,8,2,0:this%jmax), sumS(nf,8,2,0:this%jmax)
      complex(kind=dbl), intent(inout) :: cr(nf,*)
    end subroutine forward_legesum_8_sub
    
    module pure subroutine mmset_16_sub(this, m, sinx, pmm, pmj2, pmj1, pmj0)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: m
      real(kind=dbl), intent(in)    :: sinx(16)
      real(kind=dbl), intent(inout) :: pmm(16)
      real(kind=dbl), intent(out)   :: pmj2(16), pmj1(16), pmj0(16)
    end subroutine mmset_16_sub
    
    module pure subroutine recursion_16_sub(this, mj, cosx, pmj2, pmj1, pmj0)
      class(T_legep), intent(in)    :: this
      integer,        intent(in)    :: mj
      real(kind=dbl), intent(in)    :: cosx(16)
      real(kind=dbl), intent(inout) :: pmj2(16), pmj1(16), pmj0(16)
    end subroutine recursion_16_sub
    
    module pure subroutine backward_sum_16_sub(this, nb, legep, cc, legesum)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nb
      real(kind=dbl),    intent(in)    :: legep(16)
      complex(kind=dbl), intent(in)    :: cc(nb)
      complex(kind=dbl), intent(inout) :: legesum(16,nb)
    end subroutine backward_sum_16_sub
    
    module pure subroutine backward_rcb_16_sub(this, nb, sumsym, sumasym, sumN, sumS)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nb
      complex(kind=dbl), intent(in)  :: sumsym(16,nb), sumasym(16,nb)
      real(kind=dbl),    intent(out) :: sumN(nb,16,2), sumS(nb,16,2)
    end subroutine backward_rcb_16_sub
    
    module pure subroutine forward_sum_16_sub(this, nf, legep, legesum, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: nf
      real(kind=dbl),    intent(in)    :: legep(16)
      complex(kind=dbl), intent(in)    :: legesum(16,nf)
      complex(kind=dbl), intent(inout) :: cr(nf)
    end subroutine forward_sum_16_sub
    
    module pure subroutine forward_rcb_16_sub(this, nf, w, sumN, sumS, sumsym, sumasym)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: nf
      real(kind=dbl),    intent(in)  :: w(16)
      real(kind=dbl),    intent(in)  :: sumN(nf,16,2), sumS(nf,16,2)
      complex(kind=dbl), intent(out) :: sumsym(16,nf), sumasym(16,nf)
    end subroutine forward_rcb_16_sub
    
    module pure subroutine backward_legesum_16_sub(this, it, nb, cc, sumN, sumS)
      class(T_legep),    intent(in)  :: this
      integer,           intent(in)  :: it, nb
      complex(kind=dbl), intent(in)  :: cc(nb,*)
      real(kind=dbl),    intent(out) :: sumN(nb,16,2,0:this%jmax), sumS(nb,16,2,0:this%jmax)
    end subroutine backward_legesum_16_sub
    
    module pure subroutine forward_legesum_16_sub(this, it, nf, sumN, sumS, cr)
      class(T_legep),    intent(in)    :: this
      integer,           intent(in)    :: it, nf
      real(kind=dbl),    intent(in)    :: sumN(nf,16,2,0:this%jmax), sumS(nf,16,2,0:this%jmax)
      complex(kind=dbl), intent(inout) :: cr(nf,*)
    end subroutine forward_legesum_16_sub
  end interface
  
end module lege_poly