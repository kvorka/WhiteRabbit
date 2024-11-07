submodule (lege_poly) init
  implicit none; contains
  
  module procedure init_lege_sub
    
    this%nLege = nLege; call this%roots_sub()
    this%jmax  = jmax; call this%coeffs_sub()
    
    this%rw(:,3) = this%rw(:,3) / wfac
    
  end procedure init_lege_sub
  
  module procedure deallocate_lege_sub
    
    if ( allocated(this%rw) ) deallocate( this%rw )
    if ( allocated(this%ab) ) deallocate( this%ab )
    
  end procedure deallocate_lege_sub
  
end submodule init