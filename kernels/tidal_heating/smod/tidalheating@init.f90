submodule (tidalheating) init
  implicit none; contains
  
  module procedure init_tidalHeating_sub
    
    this%nd = nd
    this%jms = jms
    
    allocate( this%htide(this%nd,this%jms) )
      this%htide = czero
    
  end procedure init_tidalHeating_sub
  
  module procedure deallocate_tidalHeating_sub
    
    if ( allocated(this%htide) ) deallocate( this%htide )
    
  end procedure deallocate_tidalHeating_sub
  
end submodule init