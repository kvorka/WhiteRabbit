submodule (gravity) init
  implicit none; contains
  
  module procedure init_gravity_sub
    
    this%gmod = gmod
    this%g    = g
    
  end procedure init_gravity_sub
  
  module procedure set_gravity_sub
    
    if ( present(Dcrust) ) this%Dcrust = Dcrust
    if ( present(omega)  ) this%omega  = omega
    if ( present(exc)    ) this%exc    = exc
        
  end procedure set_gravity_sub
  
  module procedure deallocate_gravity_sub
  end procedure deallocate_gravity_sub
  
end submodule init