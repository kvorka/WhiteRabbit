submodule (IceCrustMod) Init_iceCrust
  implicit none; contains
  
  module subroutine init_iceCrust_sub(this)
    class(T_iceCrust), intent(inout) :: this
   
    call this%init_ice_sub(jmax_in = jmax_ice, rheol_in = 'viscel', n_iter = n_iter_ice)
    
    this%cf = one
    
    call this%init_eq_temp_sub( rhs=.true. , nl=.true.  )
    call this%init_eq_mech_sub( rhs=.true. , nl=.false. )
    
    call this%tdheat%init_sub(this%nd, this%jms)
    call this%bnd%init_temp_up_sub()
      call this%set_surfTemp_sub()
    
    call this%mparams%init_visc_radial_sub()
    call this%mparams%init_cp_radial_sub()
    call this%mparams%init_lambda_radial_sub()
    call this%mparams%init_alpha_radial_sub()
    
    call this%tides%init_sub()
    
  end subroutine init_iceCrust_sub
  
end submodule Init_iceCrust