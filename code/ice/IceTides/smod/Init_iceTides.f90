submodule (IceTidesMod) Init_iceTides
  implicit none; contains
  
  module subroutine init_iceTides_sub(this)
    class(T_iceTides), intent(inout) :: this
    
    call this%init_ice_sub(jmax_in=2, rheol_in='viscel', n_iter=n_iter_tides, noharm=.true.)
    
    this%cf = one
    
    call this%init_eq_mech_sub( rhs=.true. , nl=.false. )
    
    call this%tdheat%init_sub(this%nd, jms4)
    call this%mparams%init_visc_radial_sub()
    
  end subroutine init_iceTides_sub
  
end submodule Init_iceTides