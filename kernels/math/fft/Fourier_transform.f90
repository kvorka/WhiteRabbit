module Fourier_transform
  !Author of the original code: Keiichi Ishioka
  !Original work: fxpack (ISPACK FORTRAN SUBROUTINE LIBRARY FOR SCIENTIFIC COMPUTING)
  use Math
  implicit none
  
  type, public :: T_fft
    integer                     :: n
    integer,        allocatable :: it(:)
    real(kind=dbl), allocatable :: t(:)
    
    contains
    
    procedure :: init_sub       => fft_init_sub
    procedure :: exec_r2c_sub   => fft_r2c_exec_sub
    procedure :: exec_c2r_sub   => fft_c2r_exec_sub
    procedure :: deallocate_sub => fft_deallocate_sub
    
    procedure, private :: fft_r2c_sub
    procedure, private :: fft_c2r_sub
    
  end type T_fft
  
  integer,        parameter :: imm = -2e4
  real(kind=dbl), parameter :: C31 = -0.5_dbl
  real(kind=dbl), parameter :: C32 = +0.86602540378443864676_dbl
  real(kind=dbl), parameter :: C51 = +0.25_dbl
  real(kind=dbl), parameter :: C52 = +0.5590169943749474241_dbl
  real(kind=dbl), parameter :: C53 = +0.6180339887498948482_dbl
  real(kind=dbl), parameter :: C54 = -0.9510565162951535721_dbl
  
  interface
    module pure subroutine fft_init_sub(this, n)
      class(T_fft), intent(inout) :: this
      integer,      intent(in)    :: n
    end subroutine fft_init_sub
    
    module pure subroutine fft_deallocate_sub(this)
      class(T_fft), intent(inout) :: this
    end subroutine fft_deallocate_sub
    
    module pure subroutine fft_r2c_exec_sub(this, m, np, x, cx)
      class(T_fft),      intent(in)    :: this
      integer,           intent(in)    :: m, np
      real(kind=dbl),    intent(inout) :: x(m,this%n)
      complex(kind=dbl), intent(out)   :: cx(m,np)
    end subroutine fft_r2c_exec_sub
    
    pure module subroutine fft_r2c_sub(this, m, x)
      class(T_fft),      intent(in)    :: this
      integer,           intent(in)    :: m
      real(kind=dbl),    intent(inout) :: x(m,2,0:this%n/2-1)
    end subroutine fft_r2c_sub
    
    module pure subroutine fft_c2r_exec_sub(this, m, np, cx, x)
      class(T_fft),      intent(in)  :: this
      integer,           intent(in)  :: m, np
      complex(kind=dbl), intent(in)  :: cx(m,np)
      real(kind=dbl),    intent(out) :: x(m,this%n)
    end subroutine fft_c2r_exec_sub
    
    module pure subroutine fft_c2r_sub(this, m, x)
      class(T_fft),   intent(in)    :: this
      integer,        intent(in)    :: m
      real(kind=dbl), intent(inout) :: x(m,2,0:this%n/2-1)
    end subroutine fft_c2r_sub
  end interface
  
  interface
    module pure subroutine fxzini(n, it, t)
      integer,        intent(in)  :: n
      integer,        intent(out) :: it(n)
      real(kind=dbl), intent(out) :: t(2,0:n-1)
    end subroutine fxzini
    
    module pure recursive subroutine fxztal(m, k, l, x, t, ic, itsum, is, j1, it1)
      integer,        intent(in)    :: m, k, l, ic, itsum, is, j1, it1
      real(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(m,2,0:*)
    end subroutine fxztal
    
    module pure subroutine fxzm2a(m, k, l, x, t)
      integer,        intent(in)    :: m, k, l
      real(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(m,2,l/2,0:1,0:k-1)
    end subroutine fxzm2a
    
    module pure subroutine fxzm2b(m, l, x)
      integer,        intent(in)    :: m, l
      real(kind=dbl), intent(inout) :: x(m,2,l/2,0:1)
    end subroutine fxzm2b
    
    module pure subroutine fxzm3a(m, k, l, x, t)
      integer,        intent(in)    :: m, k, l
      real(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(m,2,l/3,0:2,0:k-1)
    end subroutine fxzm3a
    
    module pure subroutine fxzm3b(m, l, x)
      integer,        intent(in)    :: m, l
      real(kind=dbl), intent(inout) :: x(m,2,l/3,0:2)
    end subroutine fxzm3b
    
    module pure subroutine fxzm4a(m, k, l, x, t)
      integer,        intent(in)    :: m, k, l
      REAL(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(m,2,l/4,0:3,0:k-1)
    end subroutine fxzm4a
    
    module pure subroutine fxzm4b(m, l, x)
      integer,        intent(in)    :: m, l
      real(kind=dbl), intent(inout) :: x(m,2,l/4,0:3)
    end subroutine fxzm4b
    
    module pure subroutine fxzm5a(m, k, l, x, t)
      integer,        intent(in)    :: m, k, l
      real(kind=dbl), intent(in)    :: t(2,0:*)
      real(kind=dbl), intent(inout) :: x(m,2,l/5,0:4,0:k-1)
    end subroutine fxzm5a
    
    module pure subroutine fxzm5b(m, l, x)
      integer,        intent(in)    :: m, l
      real(kind=dbl), intent(inout) :: x(m,2,l/5,0:4)
    end subroutine fxzm5b
  end interface
  
end module Fourier_transform