module IceViscosity
  use Math
  implicit none
  
  public  :: andrade_visc_fn, goldsby_visc_fn
  private :: vdiff_visc_fn, gbdiff_visc_fn, disl_visc_fn, basal_visc_fn, gbs_visc_fn
  
  contains
  
  pure real(kind=dbl) function andrade_visc_fn(mu, angfreq, visc)
    real(kind=dbl), parameter  :: and_a  = 0.33_dbl
    real(kind=dbl), parameter  :: cosgam = gamma(1+and_a) * cos( and_a * pi / 2 )
    real(kind=dbl), parameter  :: singam = gamma(1+and_a) * sin( and_a * pi / 2 )
    real(kind=dbl), intent(in) :: mu, angfreq, visc
    
    andrade_visc_fn = visc *  ( (1 + ( mu / ( angfreq * visc ) )**(and_a  ) * cosgam ) / &
                              & (1 + ( mu / ( angfreq * visc ) )**(and_a-1) * singam )   )
    
  end function andrade_visc_fn
  
  pure real(kind=dbl) function goldsby_visc_fn(diam, temperature, stress)
    real(kind=dbl), intent(in) :: diam, temperature, stress
    
    if ( stress > 0._dbl) then
      goldsby_visc_fn = 1 / ( 1 / vdiff_visc_fn(diam, temperature)        +                                       &
                            & 1 / gbdiff_visc_fn(diam, temperature)       +                                       &
                            & 1 / disl_visc_fn(diam, temperature, stress) +                                       &
                            & 1 / ( basal_visc_fn(temperature, stress) + gbs_visc_fn(diam, temperature, stress) ) )
    else
      goldsby_visc_fn = 1 / ( 1 / vdiff_visc_fn(diam, temperature) + &
                            & 1 / gbdiff_visc_fn(diam, temperature)  )
    end if                            
    
  end function goldsby_visc_fn
  
  pure real(kind=dbl) function goldsby_diffvisc_fn(diam, temperature)
    real(kind=dbl), intent(in) :: diam, temperature
    
    goldsby_diffvisc_fn = 1 / ( 1 / vdiff_visc_fn(diam, temperature) + 1 / gbdiff_visc_fn(diam, temperature) )
    
  end function goldsby_diffvisc_fn
  
    pure real(kind=dbl) function vdiff_visc_fn(diam, temperature)
      real(kind=dbl), intent(in) :: diam, temperature
      real(kind=dbl)             :: a, e
      
      a = 2 * 1.4d-7
      e = 5.9d4
      
      vdiff_visc_fn = ( temperature * diam**2 / a ) * exp( e / ( rgas * temperature ) )
      
    end function vdiff_visc_fn
    
    pure real(kind=dbl) function gbdiff_visc_fn(diam, temperature)
      real(kind=dbl), intent(in) :: diam, temperature
      real(kind=dbl)             :: a, e
      
      a = 2 * 3.0d-16
      e = 4.9d4
      
      gbdiff_visc_fn = ( temperature * diam**3 / a ) * exp( e / ( rgas * temperature ) )
      
    end function gbdiff_visc_fn
    
    pure real(kind=dbl) function disl_visc_fn(diam, temperature, stress)
      real(kind=dbl), intent(in) :: diam, temperature, stress
      real(kind=dbl)             :: a, e
      
      if (temperature > 258) then
        a = 2 * 4.7d5
        e = 1.8d5
      else
        a = 2 * 3.1d-18
        e = 6.0d4
      end if
      
      disl_visc_fn = exp( e / ( rgas * temperature ) ) / ( a * stress**3 )
      
    end function disl_visc_fn
    
    pure real(kind=dbl) function basal_visc_fn(temperature, stress)
      real(kind=dbl), intent(in) :: temperature, stress
      real(kind=dbl)             :: a, e
      
      a = 2 * 7.1d-7
      e = 6.0d4
      
      basal_visc_fn = exp( e / ( rgas * temperature ) ) / ( a * stress**1.4_dbl )
      
    end function basal_visc_fn
    
    pure real(kind=dbl) function gbs_visc_fn(diam, temperature, stress)
      real(kind=dbl), intent(in) :: diam, temperature, stress
      real(kind=dbl)             :: a, e
      
      if ( temperature > 255 ) then
        a = 2 * 1.1d16
        e = 1.92d5
      else
        a = 2 * 1.4d-13
        e = 4.9d4
      end if
      
      gbs_visc_fn = ( diam**1.4_dbl / ( a * stress**0.8_dbl ) ) * exp( e / ( rgas * temperature ) )
      
    end function gbs_visc_fn
  
end module IceViscosity