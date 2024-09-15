submodule(IceMod) ParametersIce
  implicit none ; contains
  
  module pure real(kind=dbl) function lambda_ice_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    real(kind=dbl)            :: lambdaI
    
    if ( this%rad_grid%r(i) < this%ru - this%hC ) then
      lambdaI = name_conductivity_fn( this%temperature_ice_r_fn(i) )
    else
      lambdaI = this%lambdaC
    end if
    
    lambda_ice_fn = lambdaI / this%lambdaU
    
  end function lambda_ice_fn
  
  module pure real(kind=dbl) function cp_ice_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    
    cp_ice_fn = name_capacity_fn( this%temperature_ice_rr_fn(i) ) / this%cU
    
  end function cp_ice_fn
  
  module pure real(kind=dbl) function alpha_ice_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    
    alpha_ice_fn = name_expansivity_fn( this%temperature_ice_rr_fn(i) ) / this%alphaU
    
  end function alpha_ice_fn
  
  module pure real(kind=dbl) function visc_ice_fn(this, i)
    class(T_ice),      intent(in)  :: this
    integer,           intent(in)  :: i
    real(kind=dbl)                 :: visc
    
    if ( this%mparams%initvisc ) then
      visc_ice_fn = c2r_fn( this%mparams%visc(1,i) ) / s4pi
    
    else 
      visc = min( goldsby_visc_fn( this%diam, this%temperature_ice_r_fn(i), this%devstress_ice_r_fn(i) ), this%cutoff )
      
      if ( this%andrade ) then
        visc_ice_fn = andrade_visc_fn(this%mu, this%omega, visc) / this%viscU
      else
        visc_ice_fn = visc / this%viscU
      end if
    end if
    
  end function visc_ice_fn
  
end submodule ParametersIce