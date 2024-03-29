submodule(IceMod) ParametersIce
  implicit none ; contains
  
  module pure real(kind=dbl) function lambda_ice_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    real(kind=dbl)            :: temp, lambdaI

    if ( this%rad_grid%r(i) < this%rad_grid%r(this%nd) - this%hC ) then
       temp = this%Tu + ( this%Td-this%Tu ) * c2r_fn( this%rad_grid%c(i,-1) * this%sol%temp_fn(i  ,1) + &
                                                    & this%rad_grid%c(i,+1) * this%sol%temp_fn(i+1,1)   ) / s4pi
       
       lambdaI = 0.4685_dbl + 488.12_dbl / temp
    else
       lambdaI = this%lambdaC
    end if
     
    lambda_ice_fn = lambdaI / this%lambdaU
    
  end function lambda_ice_fn
  
  module pure real(kind=dbl) function cp_ice_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    
    cp_ice_fn = ( 185._dbl + 7.037_dbl * ( this%Tu + (this%Td-this%Tu) * c2r_fn( this%sol%temp_fn(i,1) ) / s4pi ) ) / this%cU
    
  end function cp_ice_fn
  
  module pure real(kind=dbl) function alpha_ice_fn(this, i)
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
    real(kind=dbl)            :: temp, v, dv
    
    temp = this%Tu + (this%Td-this%Tu) * c2r_fn( this%sol%temp_fn(i,1) ) / sqrt(4*pi)
    
    v  = A0 +    A3 * temp**3 +     A4 * temp**4 +     A5 * temp**5 +     A6 * temp**6 +     A7 * temp**7
    dv =      3* A3 * temp**2 + 4 * A4 * temp**3 + 5 * A5 * temp**4 + 6 * A6 * temp**5 + 7 * A7 * temp**6
    
    alpha_ice_fn = ( dv / v ) / this%alphaU
    
  end function alpha_ice_fn
  
  module pure real(kind=dbl) function visc_ice_fn(this, i)
    class(T_ice),      intent(in)  :: this
    integer,           intent(in)  :: i
    real(kind=dbl)                 :: fac1, fac2, factemp, facstress, temp, stress, visc
    complex(kind=dbl), allocatable :: devtens(:)
    
    fac1 = this%rad_grid%c(i,-1)
    fac2 = this%rad_grid%c(i,+1)
    
    factemp   = this%Td-this%Tu
    facstress = this%viscU * this%kappaU / this%D_ud**2
    
    allocate( devtens(this%jmt) ) ; devtens = this%sol%deviatoric_stress_jml2_fn(i)
      
      temp   = factemp   * c2r_fn( fac1 * this%sol%temp_fn(i,1) + fac2 * this%sol%temp_fn(i+1,1) ) / s4pi + this%Tu
      stress = facstress * sqrt( tensproduct_fn(this%jmax, devtens, devtens) ) / s4pi
      
    deallocate( devtens )
    
    visc = min( goldsby_visc_fn(this%diam, temp, stress), this%cutoff )
      if ( .not. this%andrade ) then
        visc_ice_fn = visc / this%viscU
      else
        visc_ice_fn = andrade_visc_fn(this%mu, this%omega, visc) / this%viscU
      end if
    
  end function visc_ice_fn
  
end submodule ParametersIce