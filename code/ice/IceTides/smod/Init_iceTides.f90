submodule (IceTidesMod) Init_iceTides
  implicit none; contains
  
  module subroutine init_iceTides_sub(this, latvisc)
    class(T_iceTides), intent(inout) :: this
    logical,           intent(in)    :: latvisc
    
    call this%init_ice_sub(jmax_in=2, rheol_in='viscel', n_iter=n_iter_tides, noharm=.true.)
      this%cf = one
    
    call this%init_eq_mech_sub( rhs=.true. , nl=.false. )
    
    if ( latvisc ) then
      call this%mparams%init_visc_sub()
    else
      call this%mparams%init_visc_radial_sub()
    end if
    
  end subroutine init_iceTides_sub
  
end submodule Init_iceTides