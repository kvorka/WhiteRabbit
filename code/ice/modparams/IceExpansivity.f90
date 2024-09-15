module IceExpansivity
  use Math
  implicit none
  
  real(kind=dbl), parameter, private :: A0 = 128.2147d0
  real(kind=dbl), parameter, private :: A3 = -1.3152d-6
  real(kind=dbl), parameter, private :: A4 = +2.4837d-8
  real(kind=dbl), parameter, private :: A5 = -1.6064d-10
  real(kind=dbl), parameter, private :: A6 = +4.6097d-13
  real(kind=dbl), parameter, private :: A7 = -4.9661d-16
  
  public :: name_expansivity_fn
  
  contains
  
  elemental pure real(kind=dbl) function name_expansivity_fn(temperature)
    real(kind=dbl), intent(in) :: temperature
    real(kind=dbl)             :: v, dv
    
    v = A0 + A3 * temperature**3 + &
      &      A4 * temperature**4 + &
      &      A5 * temperature**5 + &
      &      A6 * temperature**6 + &
      &      A7 * temperature**7
    
    dv = 3 * A3 * temperature**2 + &
       & 4 * A4 * temperature**3 + &
       & 5 * A5 * temperature**4 + &
       & 6 * A6 * temperature**5 + &
       & 7 * A7 * temperature**6
    
    name_expansivity_fn = dv / v
    
  end function name_expansivity_fn
  
end module IceExpansivity