submodule (icetides) init
  implicit none; contains
  
  module procedure init_iceTides_sub
    
    call this%init_ice_sub(jmax_in=2, rheol_in='viscel', n_iter=n_iter_tides, noharm=.true.)
    
    this%cf = one
    
    call this%init_eq_mech_sub()
    
    call this%tdheat%init_sub(this%nd, jms4)
    call this%mparams%init_visc_radial_sub()
    
  end procedure init_iceTides_sub
  
end submodule init