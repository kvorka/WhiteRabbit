submodule (lege_poly) init
  implicit none; contains
  
  module procedure init_lege_sub
    
    this%nLege = nLege; call this%roots_sub()
    this%jmax  = jmax; call this%coeffs_sub()
    
    this%rootsweights(3,:) = this%rootsweights(3,:) / wfac
    
  end procedure init_lege_sub
  
  module procedure deallocate_lege_sub
    
    if ( allocated(this%rootsweights) ) deallocate( this%rootsweights )
    
    if ( allocated(this%amj) ) deallocate( this%amj )
    if ( allocated(this%bmj) ) deallocate( this%bmj )
    if ( allocated(this%cmm) ) deallocate( this%cmm )
    
  end procedure deallocate_lege_sub
  
end submodule init