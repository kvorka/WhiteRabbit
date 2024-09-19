submodule(IceMod) ParametersIce
  implicit none ; contains
  
  module pure real(kind=dbl) function lambda_r_ice_fn(this, ir)
    class(T_ice), intent(in) :: this
    integer,      intent(in) :: ir
    
    if ( this%rad_grid%r(ir) < this%ru - this%hC ) then
      
      if ( this%mparams%initlambda ) then
        lambda_r_ice_fn = c2r_fn( this%mparams%lambda(1,ir) ) / s4pi
      else
        lambda_r_ice_fn = name_conductivity_fn( this%average_temperature_ice_ir_fn(ir) ) / this%lambdaU
      end if
      
    else
      
      lambda_r_ice_fn = this%lambdaC / this%lambdaU
      
    end if
    
  end function lambda_r_ice_fn
  
  module pure real(kind=dbl) function lambda_rr_ice_fn(this, ir)
    class(T_ice), intent(in) :: this
    integer,      intent(in) :: ir
    
    if ( this%rad_grid%rr(ir) < this%ru - this%hC ) then
      
      if ( this%mparams%initlambda ) then
        lambda_rr_ice_fn = c2r_fn( this%rad_grid%cc(ir,-1) * this%mparams%lambda(1,ir-1) + &
                                 & this%rad_grid%cc(ir,+1) * this%mparams%lambda(1,ir  )   ) / s4pi
      else
        lambda_rr_ice_fn = name_conductivity_fn( this%average_temperature_ice_irr_fn(ir) ) / this%lambdaU
      end if
      
    else
      
      lambda_rr_ice_fn = this%lambdaC / this%lambdaU
      
    end if
    
  end function lambda_rr_ice_fn
  
  module pure real(kind=dbl) function cp_r_ice_fn(this, ir)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: ir
    
    cp_r_ice_fn = name_capacity_fn( this%average_temperature_ice_ir_fn(ir) ) / this%cU
    
  end function cp_r_ice_fn
  
  module pure real(kind=dbl) function cp_rr_ice_fn(this, ir)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: ir
    
    cp_rr_ice_fn = name_capacity_fn( this%average_temperature_ice_irr_fn(ir) ) / this%cU
    
  end function cp_rr_ice_fn
  
  module pure real(kind=dbl) function alpha_r_ice_fn(this, ir)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: ir
    
    alpha_r_ice_fn = name_expansivity_fn( this%average_temperature_ice_ir_fn(ir) ) / this%alphaU
    
  end function alpha_r_ice_fn
  
  module pure real(kind=dbl) function alpha_rr_ice_fn(this, ir)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: ir
    
    alpha_rr_ice_fn = name_expansivity_fn( this%average_temperature_ice_irr_fn(ir) ) / this%alphaU
    
  end function alpha_rr_ice_fn
  
  module pure real(kind=dbl) function visc_r_ice_fn(this, ir)
    class(T_ice),      intent(in)  :: this
    integer,           intent(in)  :: ir
    real(kind=dbl)                 :: visc
    
    if ( this%mparams%initvisc ) then
      visc_r_ice_fn = c2r_fn( this%mparams%visc(1,ir) ) / s4pi / this%viscU
    
    else 
      visc = min( goldsby_visc_fn( this%diam, this%average_temperature_ice_ir_fn(ir), &
                                 & this%average_stress_ice_ir_fn(ir) ), this%cutoff   )
      
      if ( this%andrade ) then
        visc_r_ice_fn = andrade_visc_fn(this%mu, this%omega, visc) / this%viscU
      else
        visc_r_ice_fn = visc / this%viscU
      end if
    end if
    
  end function visc_r_ice_fn
  
  module pure real(kind=dbl) function visc_rr_ice_fn(this, ir)
    class(T_ice),      intent(in)  :: this
    integer,           intent(in)  :: ir
    real(kind=dbl)                 :: visc
    
    if ( this%mparams%initvisc ) then
      visc_rr_ice_fn = c2r_fn( this%rad_grid%cc(ir,-1) * this%mparams%visc(1,ir-1) + &
                             & this%rad_grid%cc(ir,+1) * this%mparams%visc(1,ir  ) ) / s4pi
    
    else 
      visc = min( goldsby_visc_fn( this%diam, this%average_temperature_ice_irr_fn(ir), this%average_stress_ice_irr_fn(ir) ), this%cutoff )
      
      if ( this%andrade ) then
        visc_rr_ice_fn = andrade_visc_fn(this%mu, this%omega, visc) / this%viscU
      else
        visc_rr_ice_fn = visc / this%viscU
      end if
    end if
    
  end function visc_rr_ice_fn
  
end submodule ParametersIce