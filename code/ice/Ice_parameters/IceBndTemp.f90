module IceBndTemp
  use math
  implicit none; public; contains
  
  pure real(kind=dbl) function name_surfaceTemp_fn(theta)
    real(kind=dbl), intent(in) :: theta
    
    name_surfaceTemp_fn = 90._dbl * ( sin(theta) )**(0.25_dbl) + 30._dbl
    
  end function name_surfaceTemp_fn
  
  pure real(kind=dbl) function name_meltingTemp_fn( thickness )
    real(kind=dbl), intent(in) :: thickness
    
    name_meltingTemp_fn = 557.2_dbl - 273 * exp( (300 + 1.242 * thickness/1e3 )**2 / 2270000 )
    
  end function name_meltingTemp_fn
  
end module IceBndTemp
