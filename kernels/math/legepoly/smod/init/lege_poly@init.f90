submodule (lege_poly) init
  implicit none; contains
  
  module subroutine init_lege_sub(this, jmax, nLege)
    class(T_legep), intent(inout) :: this
    integer,        intent(in)    :: jmax, nLege
    
    this%nLege = nLege; call this%roots_sub()
    this%jmax  = jmax; call this%coeffs_sub()
    
  end subroutine init_lege_sub
  
  module pure subroutine deallocate_lege_sub(this)
    class(T_legep), intent(inout) :: this
    
    if ( allocated(this%rootsweights) ) deallocate( this%rootsweights )
    
    if ( allocated(this%amj) ) deallocate( this%amj )
    if ( allocated(this%bmj) ) deallocate( this%bmj )
    if ( allocated(this%cmm) ) deallocate( this%cmm )
    
  end subroutine deallocate_lege_sub
  
end submodule init