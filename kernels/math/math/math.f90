module math
  use iso_fortran_env, only: real64, real128
  use omp_lib
  use iso_c_binding
  use alignement
  implicit none; public
  
  integer,           parameter :: dbl  = real64    !double precision
  integer,           parameter :: qbl  = real128   !quadruple precision
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
    
    module subroutine alloc_aligned1d_sub( n, c_arr, f_arr )
      integer,                 intent(in)  :: n
      type(c_ptr),             intent(out) :: c_arr
      real(kind=dbl), pointer, intent(out) :: f_arr(:)
    end subroutine alloc_aligned1d_sub
    
    module subroutine free_aligned1d_sub( c_arr, f_arr )
      type(c_ptr),             intent(inout) :: c_arr
      real(kind=dbl), pointer, intent(inout) :: f_arr(:)
    end subroutine free_aligned1d_sub
    
    module subroutine alloc_aligned2d_sub( n1, n2, c_arr, f_arr )
      integer,                 intent(in)  :: n1, n2
      type(c_ptr),             intent(out) :: c_arr
      real(kind=dbl), pointer, intent(out) :: f_arr(:,:)
    end subroutine alloc_aligned2d_sub
    
    module subroutine free_aligned2d_sub( c_arr, f_arr )
      type(c_ptr),             intent(inout) :: c_arr
      real(kind=dbl), pointer, intent(inout) :: f_arr(:,:)
    end subroutine free_aligned2d_sub
    
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
    type(c_ptr) function malloc(alignement, n) bind(C, name='aligned_alloc')
      import         :: c_ptr
      integer, value :: alignement, n
    end function malloc
    
    subroutine free(ptr) bind(C, name="free")
      import             :: c_ptr
      type(c_ptr), value :: ptr
    end subroutine free
  end interface
  
end module math
