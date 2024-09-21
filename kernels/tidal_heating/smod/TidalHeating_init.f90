submodule (Tidal_heating) TidalHeating_init
  implicit none; contains
  
  module pure subroutine init_tidalHeating_sub(this, nd, jms)
    class(T_tidalHeating), intent(inout) :: this
    integer,               intent(in)    :: nd, jms
    
    this%nd = nd
    this%jms = jms
    
    allocate( this%htide(this%nd,this%jms) )
      this%htide = czero
    
  end subroutine init_tidalHeating_sub
  
  module pure subroutine deallocate_tidalHeating_sub(this)
    class(T_tidalHeating), intent(inout) :: this
    
    if ( allocated(this%htide) ) deallocate( this%htide )
    
  end subroutine deallocate_tidalHeating_sub
  
end submodule TidalHeating_init