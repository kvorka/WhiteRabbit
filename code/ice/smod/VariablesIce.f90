submodule(IceMod) VariablesIce
  implicit none; contains
  
  module pure real(kind=dbl) function temperature_ice_r_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    
    temperature_ice_r_fn = this%Tu + ( this%Td-this%Tu ) * c2r_fn( this%rad_grid%c(i,-1) * this%sol%temp_fn(i  ,1) + &
                                                                 & this%rad_grid%c(i,+1) * this%sol%temp_fn(i+1,1)   ) / s4pi
    
  end function temperature_ice_r_fn
  
  module pure real(kind=dbl) function temperature_ice_rr_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    
    temperature_ice_rr_fn = this%Tu + ( this%Td-this%Tu ) * c2r_fn( this%sol%temp_fn(i,1) ) / s4pi
    
  end function temperature_ice_rr_fn
  
  module pure real(kind=dbl) function devstress_ice_r_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    
    devstress_ice_r_fn = this%viscU * this%kappaU / this%D_ud**2 * &
                       & sqrt( tensnorm2_fn(this%jmax, this%sol%deviatoric_stress_jml2_fn(i)) ) / s4pi
    
  end function devstress_ice_r_fn
  
end submodule VariablesIce