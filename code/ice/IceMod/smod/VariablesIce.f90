submodule (IceMod) VariablesIce
  implicit none; contains
  
  module pure real(kind=dbl) function avrg_temperature_ice_ir_fn(this, ir)
    class(T_ice), intent(in) :: this
    integer,      intent(in) :: ir
    
    avrg_temperature_ice_ir_fn = this%Tu + ( this%Td-this%Tu ) * c2r_fn( this%temp_r_fn(ir,1) ) / s4pi
    
  end function avrg_temperature_ice_ir_fn
  
  module pure real(kind=dbl) function avrg_temperature_ice_irr_fn(this, ir)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: ir
    
    avrg_temperature_ice_irr_fn = this%Tu + ( this%Td-this%Tu ) * c2r_fn( this%temp_rr_fn(ir,1) ) / s4pi
    
  end function avrg_temperature_ice_irr_fn
  
  module pure real(kind=dbl) function avrg_stress_ice_ir_fn(this, ir)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: ir
    real(kind=dbl)            :: fac
    
    fac = this%viscU * this%kappaU / this%D_ud**2
    
    avrg_stress_ice_ir_fn = fac * sqrt( tensnorm2_fn(this%jmax, this%sol%deviatoric_stress_jml2_fn(ir)) ) / s4pi
    
  end function avrg_stress_ice_ir_fn
  
  module pure real(kind=dbl) function avrg_stress_ice_irr_fn(this, ir)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: ir
    
    avrg_stress_ice_irr_fn = zero
    
  end function avrg_stress_ice_irr_fn
  
end submodule VariablesIce