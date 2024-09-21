module IceSurfaceTemp
  use Math
  implicit none; public
  
  contains
  
  pure elemental real(kind=dbl) function name_surfaceTemp_fn(theta, phi)
    real(kind=dbl), intent(in) :: theta, phi
    
    !!***********************!!
    !!  Treba pridat nejaku  !!
    !!  implementaciu.       !!
    !!***********************!!
    
    name_surfaceTemp_fn = zero
    
  end function name_surfaceTemp_fn
  
end module IceSurfaceTemp
