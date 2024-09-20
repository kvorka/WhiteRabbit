module IceCapacity
  use Math
  implicit none; public
  
  real(kind=dbl), parameter, private :: c1 = 185._dbl
  real(kind=dbl), parameter, private :: c2 = 7.037_dbl
  
  public :: name_capacity_fn
  
  contains
  
  pure elemental real(kind=dbl) function name_capacity_fn(temperature)
    real(kind=dbl), intent(in) :: temperature
    
    name_capacity_fn = c1 + c2 * temperature
    
  end function name_capacity_fn
  
end module IceCapacity