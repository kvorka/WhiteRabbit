module ice_bndtemp
  use math
  use ice_constants
  implicit none; public
  
  real(kind=dbl), parameter, private :: sun_lumin = 3.830d26
  real(kind=dbl), parameter, private :: stfbtzconst = 5.670d-8
  real(kind=dbl), parameter, private :: Lconst = (1-albedo) * sun_lumin / ( 8 * pi**3 * sun_to_body_dist**2 * sqrt(1-exc**2) )
  
  contains
  
  pure real(kind=dbl) function ward_phi_integral_fn(theta)
    real(kind=dbl), intent(in) :: theta
    integer,        parameter  :: nphi = 1000
    real(kind=dbl), parameter  :: dphi = 2*pi / (nphi-1)
    integer                    :: iph
    real(kind=dbl)             :: integrand, phi
    
    ward_phi_integral_fn = 0._dbl
    
    do iph = 1, nphi-1
      phi       = (iph-half) / (nphi-1) * pi / 180
      integrand = sqrt( 1 - ( cos(theta)*cos(obl) - sin(theta)*sin(exc)*sin(phi) )**2 )
      
      ward_phi_integral_fn = ward_phi_integral_fn + integrand * dphi
    end do
    
  end function ward_phi_integral_fn
  
  pure real(kind=dbl) function name_surfaceTemp_fn(theta, q)
    real(kind=dbl), intent(in) :: theta, q
    
    name_surfaceTemp_fn = ( ( Lconst * ward_phi_integral_fn(theta) - q ) / stfbtzconst )**(0.25_dbl)
    
  end function name_surfaceTemp_fn
  
  pure real(kind=dbl) function name_meltingTemp_fn( thickness )
    real(kind=dbl), intent(in) :: thickness
    
    name_meltingTemp_fn = 557.2_dbl - 273 * exp( (300 + 1.242 * thickness/1e3 )**2 / 2270000 )
    
  end function name_meltingTemp_fn
  
end module ice_bndtemp
