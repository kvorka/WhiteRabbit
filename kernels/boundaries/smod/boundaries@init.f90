submodule (boundaries) init
  implicit none; contains
  
  module procedure init_boundaries_sub
    
    this%jms = jms
    
  end procedure init_boundaries_sub
  
  module procedure deallocate_boundaries_sub
    
    if ( allocated(this%temp_up) ) deallocate( this%temp_up )
    if ( allocated(this%flux_up) ) deallocate( this%flux_up )
    
    if ( allocated(this%u_dn) ) deallocate( this%u_dn )
    if ( allocated(this%u_up) ) deallocate( this%u_up )
    if ( allocated(this%u_I2) ) deallocate( this%u_I2 )
    if ( allocated(this%u_C ) ) deallocate( this%u_C  )
    
    if ( allocated(this%t_dn) ) deallocate( this%t_dn )
    if ( allocated(this%t_up) ) deallocate( this%t_up )
    
    if ( allocated(this%v_up) ) deallocate( this%v_up )
    if ( allocated(this%v_dn) ) deallocate( this%v_dn )
    
  end procedure deallocate_boundaries_sub
  
end submodule init