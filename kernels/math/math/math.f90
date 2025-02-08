module math
  use iso_fortran_env, only: real64, real128
  use omp_lib
  use iso_c_binding
  implicit none; public
  
  integer,           parameter :: dbl   = real64
  integer,           parameter :: qbl   = real128
  integer,           parameter :: step  = 16
  real(kind=dbl),    parameter :: deps  = 1.0d-15
  real(kind=qbl),    parameter :: qeps  = 1.0d-28
  real(kind=dbl),    parameter :: zero  = 0._dbl
  real(kind=qbl),    parameter :: qzero = 0._qbl
  real(kind=dbl),    parameter :: half  = 0.5_dbl
  real(kind=dbl),    parameter :: one   = 1._dbl
  real(kind=dbl),    parameter :: two   = 2._dbl
  real(kind=qbl),    parameter :: qone  = 1._qbl
  real(kind=dbl),    parameter :: sq2_1 = 1 / sqrt(2._dbl)
  real(kind=dbl),    parameter :: pi    = acos(-one)
  real(kind=qbl),    parameter :: qpi   = acos(-qone)
  real(kind=dbl),    parameter :: s4pi  = sqrt(4*pi)
  complex(kind=dbl), parameter :: cone  = cmplx(one,  zero, kind=dbl)
  complex(kind=dbl), parameter :: cunit = cmplx(zero, one , kind=dbl)
  complex(kind=dbl), parameter :: czero = cmplx(zero, zero, kind=dbl)
  complex(kind=dbl), parameter :: cs4pi = cmplx(s4pi, zero, kind=dbl)
  
  interface
    module elemental function int2str_fn(n) result(str)
      integer,          intent(in) :: n
      character(len=10)            :: str
    end function int2str_fn
    
    module elemental real(kind=dbl) function i2r_fn(ix)
      integer, intent(in) :: ix
    end function i2r_fn
    
    module elemental complex(kind=dbl) function r2c_fn(x)
      real(kind=dbl), intent(in) :: x
    end function r2c_fn
    
    module elemental real(kind=dbl) function c2r_fn(cx)
      complex(kind=dbl), intent(in) :: cx
    end function c2r_fn
    
    module pure subroutine zero_rarray_sub(length, arr)
      integer,        intent(in)  :: length
      real(kind=dbl), intent(out) :: arr(*)
    end subroutine zero_rarray_sub
    
    module pure subroutine zero_carray_sub(length, arr)
      integer,           intent(in)  :: length
      complex(kind=dbl), intent(out) :: arr(*)
    end subroutine zero_carray_sub
  end interface
  
  interface
     subroutine calloc_sub(n1, n2, cptr) bind(C, name='memalloc')
      import                     :: c_ptr
      type(c_ptr), intent(inout) :: cptr
      integer,     intent(in)    :: n1, n2
    end subroutine calloc_sub

    subroutine cfree_sub(cptr) bind(C, name='memfree')
      import                     :: c_ptr
      type(c_ptr), intent(inout) :: cptr
    end subroutine cfree_sub
    
    subroutine czero_sub(n, arr) bind(C, name='memzero')
      import                      :: dbl
      integer,        intent(in)  :: n
      real(kind=dbl), intent(out) :: arr(*)
    end subroutine czero_sub
  end interface
  
end module math
