module Math
  use iso_fortran_env, only: real64
  implicit none
  
  integer,           parameter, public :: dbl   = real64
  real(kind=dbl),    parameter, public :: zero  = 0._dbl
  real(kind=dbl),    parameter, public :: one   = 1._dbl
  real(kind=dbl),    parameter, public :: pi    = acos(-one)
  real(kind=dbl),    parameter, public :: s4pi  = sqrt(4*pi)
  complex(kind=dbl), parameter, public :: cunit = cmplx(zero, one , kind=dbl)
  complex(kind=dbl), parameter, public :: czero = cmplx(zero, zero, kind=dbl)
  complex(kind=dbl), parameter, public :: cone  = cmplx(one , zero, kind=dbl)
  complex(kind=dbl), parameter, public :: cs4pi = cmplx(s4pi, zero, kind=dbl)
  real(kind=dbl),    parameter, public :: kappa = 6.670d-11
  real(kind=dbl),    parameter, public :: rgas  = 8.31_dbl

  public :: int2str_fn
  public :: i2r_fn
  public :: r2c_fn
  public :: c2r_fn

  contains

  function int2str_fn(n) result(str)
    integer,          intent(in) :: n
    character(len=10)            :: str 

    write(str,'(1I4)') n

  end function int2str_fn
  
  elemental real(kind=dbl) function i2r_fn(ix)
    integer, intent(in) :: ix
    
    i2r_fn = real(ix, kind=dbl)
    
  end function i2r_fn
  
  elemental complex(kind=dbl) function r2c_fn(x)
    real(kind=dbl), intent(in) :: x
    
    r2c_fn = cmplx(x, zero, kind=dbl)
    
  end function r2c_fn
  
  elemental real(kind=dbl) function c2r_fn(cx)
    complex(kind=dbl), intent(in) :: cx
    
    c2r_fn = real(cx, kind=dbl)
    
  end function c2r_fn
  
end module Math
