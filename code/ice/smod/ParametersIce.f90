submodule(IceMod) ParametersIce
  implicit none

  contains
  
  pure real(kind=dbl) function lambda_ice_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    real(kind=dbl)            :: temp, lambdaI

    if ( this%rad_grid%r(i) < this%rad_grid%r(this%nd) - this%hC ) then
       temp = this%Tu + (this%Td-this%Tu) * real( this%rad_grid%c(i,-1) * this%sol%temp_fn(i  ,1) + &
                                                & this%rad_grid%c(i,+1) * this%sol%temp_fn(i+1,1)    , kind=dbl) / sqrt(4*pi)
       
       lambdaI = 0.4685_dbl + 488.12_dbl / temp
       
    else
       lambdaI = this%lambdaC
       
    end if
     
    lambda_ice_fn = lambdaI / this%lambdaU
    
  end function lambda_ice_fn
  
  pure real(kind=dbl) function cp_ice_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    real(kind=dbl)            :: temp, cI
    
    temp = this%Tu + (this%Td-this%Tu) * real( this%sol%temp_fn(i,1), kind=dbl ) / sqrt(4*pi)
    cI = 185._dbl + 7.037_dbl * temp
    
    cp_ice_fn = cI / this%cU
    
  end function cp_ice_fn
  
  pure real(kind=dbl) function alpha_ice_fn(this, i)
    real(kind=dbl), parameter :: A0 = 128.2147_dbl
    real(kind=dbl), parameter :: A1 = 0._dbl
    real(kind=dbl), parameter :: A2 = 0._dbl
    real(kind=dbl), parameter :: A3 = -1.3152d-6
    real(kind=dbl), parameter :: A4 = +2.4837d-8
    real(kind=dbl), parameter :: A5 = -1.6064d-10
    real(kind=dbl), parameter :: A6 = +4.6097d-13
    real(kind=dbl), parameter :: A7 = -4.9661d-16
    
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    real(kind=dbl)            :: temp, v, dv, alphaI
    
    temp = this%Tu + (this%Td-this%Tu) * real( this%sol%temp_fn(i,1), kind=dbl ) / sqrt(4*pi)
    
    v  = A0 +    A3 * temp**3 +     A4 * temp**4 +     A5 * temp**5 +     A6 * temp**6 +     A7 * temp**7
    dv =      3* A3 * temp**2 + 4 * A4 * temp**3 + 5 * A5 * temp**4 + 6 * A6 * temp**5 + 7 * A7 * temp**6
    
    alphaI = dv/v
    
    alpha_ice_fn = alphaI / this%alphaU
    
  end function alpha_ice_fn
  
  pure real(kind=dbl) function visc_ice_fn(this, i)
    real(kind=dbl), parameter :: a = 9.0d-8
    real(kind=dbl), parameter :: e = 59.0d+3
    real(kind=dbl), parameter :: and_a  = 0.33_dbl
    real(kind=dbl), parameter :: cosgam = gamma(1+and_a) * cos( and_a * pi / 2 )
    real(kind=dbl), parameter :: singam = gamma(1+and_a) * sin( and_a * pi / 2 )
    
    class(T_ice),      intent(in) :: this
    integer,           intent(in) :: i
    real(kind=dbl)                :: temp, viscI
    
    temp = this%Tu + (this%Td-this%Tu) * real( this%rad_grid%c(i,-1) * this%sol%temp_fn(i,1) + &
                                             & this%rad_grid%c(i,+1) * this%sol%temp_fn(i+1,1), kind=dbl) / sqrt(4*pi)
    
    viscI = min( temp * this%diam**2 * exp(e / rgas / temp) / a / 2, this%cutoff )
    
    if ( .not. this%andrade ) then
      visc_ice_fn = viscI / this%viscU
    else
      visc_ice_fn = viscI / this%viscU * ( (1 + ( this%mu / ( this%omega * viscI ) )**(and_a  ) * cosgam ) / &
                                         & (1 + ( this%mu / ( this%omega * viscI ) )**(and_a-1) * singam )   )
    end if
    
  end function visc_ice_fn
  
end submodule ParametersIce