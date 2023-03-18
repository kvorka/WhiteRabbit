module Math
  use iso_fortran_env, only: real64
  implicit none; private
  
  integer,           parameter, public :: dbl = real64
  real(kind=dbl),    parameter, public :: pi  = acos(-1._dbl)
  complex(kind=dbl), parameter, public :: cunit = cmplx(0._dbl, 1._dbl, kind=dbl)
  complex(kind=dbl), parameter, public :: czero = cmplx(0._dbl, 0._dbl, kind=dbl)
  real(kind=dbl),    parameter, public :: kappa = 6.670d-11
  real(kind=dbl),    parameter, public :: rgas = 8.31_dbl

  public :: int2str_fn

  contains

  function int2str_fn(n) result(str)
    integer,          intent(in) :: n
    character(len=10)            :: str 

    write(str,'(1I4)') n

  end function int2str_fn
  
end module Math