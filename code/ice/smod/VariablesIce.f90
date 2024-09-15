submodule(IceMod) VariablesIce
  implicit none; contains
  
  module pure real(kind=dbl) function temperature_ice_r_fn(this, i)
    class(T_ice), intent(in) :: this
    integer,      intent(in) :: i
    
    temperature_ice_r_fn = this%Tu + ( this%Td-this%Tu ) * c2r_fn( this%temp_r_fn(i,1) ) / s4pi
    
  end function temperature_ice_r_fn
  
  module pure subroutine temperature_ice_r_jm_sub(this, i, temperature)
    class(T_ice),      intent(in)  :: this
    integer,           intent(in)  :: i
    complex(kind=dbl), intent(out) :: temperature(:)
    integer                        :: ijm
    
    do concurrent ( ijm = 1:this%jms )
      temperature(ijm) = this%temp_r_fn(i,ijm)
    end do
    
  end subroutine temperature_ice_r_jm_sub
  
  module pure real(kind=dbl) function temperature_ice_rr_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    
    temperature_ice_rr_fn = this%Tu + ( this%Td-this%Tu ) * c2r_fn( this%temp_rr_fn(i,1) ) / s4pi
    
  end function temperature_ice_rr_fn
  
  module pure real(kind=dbl) function devstress_ice_r_fn(this, i)
    class(T_ice),  intent(in) :: this
    integer,       intent(in) :: i
    
    devstress_ice_r_fn = this%viscU * this%kappaU / this%D_ud**2 * &
                       & sqrt( tensnorm2_fn(this%jmax, this%sol%deviatoric_stress_jml2_fn(i)) ) / s4pi
    
  end function devstress_ice_r_fn
  
end submodule VariablesIce