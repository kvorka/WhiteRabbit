module Math
  use iso_fortran_env, only: real64
  implicit none
  
  integer,           parameter, public :: dbl = real64
  real(kind=dbl),    parameter, public :: pi  = acos(-1._dbl)
  complex(kind=dbl), parameter, public :: cunit = cmplx(0._dbl, 1._dbl, kind=dbl)
  complex(kind=dbl), parameter, public :: czero = cmplx(0._dbl, 0._dbl, kind=dbl)
  complex(kind=dbl), parameter, public :: cone  = cmplx(1._dbl, 0._dbl, kind=dbl)
  complex(kind=dbl), parameter, public :: cs4pi = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
  real(kind=dbl),    parameter, public :: kappa = 6.670d-11
  real(kind=dbl),    parameter, public :: rgas = 8.31_dbl

  public :: int2str_fn
  public :: r2c_fn
  public :: c2r_fn

  contains

  function int2str_fn(n) result(str)
    integer,          intent(in) :: n
    character(len=10)            :: str 

    write(str,'(1I4)') n

  end function int2str_fn
  
  elemental function r2c_fn(x) result(cx)
    real(kind=dbl),   intent(in) :: x
    complex(kind=dbl)            :: cx
    
    cx = cmplx(x, 0._dbl, kind=dbl)
    
  end function r2c_fn
  
  elemental function c2r_fn(cx) result(x)
    complex(kind=dbl), intent(in) :: cx
    real(kind=dbl)                :: x
    
    x = real(cx, kind=dbl)
    
  end function c2r_fn
  
end module Math