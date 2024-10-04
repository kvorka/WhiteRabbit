module ice_cond
  use math
  implicit none; public
  
  real(kind=dbl), parameter, private :: c1 = 0.4685_dbl
  real(kind=dbl), parameter, private :: c2 = 488.12_dbl
  
  public :: name_conductivity_fn
  
  contains
  
  pure elemental real(kind=dbl) function name_conductivity_fn(temperature)
    real(kind=dbl), intent(in) :: temperature
    
    name_conductivity_fn = c1 + c2 / temperature
    
  end function name_conductivity_fn
  
end module ice_cond