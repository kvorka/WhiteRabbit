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
    
    cp_ice_fn = ( 185._dbl + 7.037_dbl * this%temperature_ice_rr_fn(i) ) / this%cU
    
  end function cp_ice_fn
  
  module pure real(kind=dbl) function alpha_ice_fn(this, i)
    real(kind=dbl), parameter :: A0 = 128.2147_dbl
    !real(kind=dbl), parameter :: A1 = zero
    !real(kind=dbl), parameter :: A2 = zero
    real(kind=dbl), parameter :: A3 = -1.3152d-6
    real(kind=dbl), parameter :: A4 = +2.4837d-8
    real(kind=dbl), parameter :: A5 = -1.6064d-10
    real(kind=dbl), parameter :: A6 = +4.6097d-13
    real(kind=dbl), parameter :: A7 = -4.9661d-16
    
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    real(kind=dbl)            :: temp, v, dv
    
    temp = this%temperature_ice_rr_fn(i)
    
    v  = A0 +    A3 * temp**3 +     A4 * temp**4 +     A5 * temp**5 +     A6 * temp**6 +     A7 * temp**7
    dv =      3* A3 * temp**2 + 4 * A4 * temp**3 + 5 * A5 * temp**4 + 6 * A6 * temp**5 + 7 * A7 * temp**6
    
    alpha_ice_fn = ( dv / v ) / this%alphaU
    
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