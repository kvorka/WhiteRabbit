submodule (ice) variables
  implicit none; contains
  
  module procedure avrg_temperature_ice_ir_fn
    
    avrg_temperature_ice_ir_fn = this%Tu + ( this%Td-this%Tu ) * c2r_fn( this%temp_r_fn(ir,1) ) / s4pi
    
  end procedure avrg_temperature_ice_ir_fn
  
  module procedure avrg_temperature_ice_irr_fn
    
    avrg_temperature_ice_irr_fn = this%Tu + ( this%Td-this%Tu ) * c2r_fn( this%temp_rr_fn(ir,1) ) / s4pi
    
  end procedure avrg_temperature_ice_irr_fn
  
  module procedure avrg_stress_ice_ir_fn
    real(kind=dbl) :: fac
    
    fac = this%viscU * this%kappaU / this%D_ud**2
    
    avrg_stress_ice_ir_fn = fac * sqrt( tensnorm2_fn(this%jmax, this%sol%deviatoric_stress_jml2_fn(ir)) ) / s4pi
    
  end procedure avrg_stress_ice_ir_fn
  
  module procedure avrg_stress_ice_irr_fn
    
    avrg_stress_ice_irr_fn = zero
    
  end procedure avrg_stress_ice_irr_fn
  
end submodule variables