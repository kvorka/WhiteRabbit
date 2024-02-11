module Math
  use iso_fortran_env, only: real64
  implicit none; public
  
  integer,           parameter :: dbl   = real64
  integer,           parameter :: step  = 16
  real(kind=dbl),    parameter :: zero  = 0._dbl
  real(kind=dbl),    parameter :: one   = 1._dbl
  real(kind=dbl),    parameter :: sq2_1 = 1 / sqrt(2._dbl)
  real(kind=dbl),    parameter :: pi    = acos(-one)
  real(kind=dbl),    parameter :: s4pi  = sqrt(4*pi)
  complex(kind=dbl), parameter :: cunit = cmplx(zero, one , kind=dbl)
  complex(kind=dbl), parameter :: czero = cmplx(zero, zero, kind=dbl)
  complex(kind=dbl), parameter :: cone  = cmplx(one , zero, kind=dbl)
  complex(kind=dbl), parameter :: cs4pi = cmplx(s4pi, zero, kind=dbl)
  real(kind=dbl),    parameter :: kappa = 6.670d-11
  real(kind=dbl),    parameter :: rgas  = 8.31_dbl
  
  integer, parameter :: addmissible_jmax(47) = [   5,   7,   9,  13,  15,  21,  27,  29,  33,  37,  45, 47,  51,  57,  61,  69,  &
                                               &  77,  87,  93,  97, 105, 117, 125, 141, 147, 157, 159, 177, 189, 197, 213, 237, &
                                               & 247, 253, 267, 285, 297, 317, 321, 357, 381, 397, 429, 447, 477, 483, 497       ]
  
end module Math
