submodule (Boundaries) Boundaries_init
  implicit none; contains
  
  module pure subroutine init_boundaries_sub(this, jms)
    class(T_boundaries), intent(inout) :: this
    integer,             intent(in)    :: jms
    
    this%jms = jms
    
  end subroutine init_boundaries_sub
  
  module pure subroutine deallocate_boundaries_sub(this)
    class(T_boundaries), intent(inout) :: this
    
    if ( allocated(this%u_dn) ) deallocate( this%u_dn )
    if ( allocated(this%u_up) ) deallocate( this%u_up )
    if ( allocated(this%u_I2) ) deallocate( this%u_I2 )
    if ( allocated(this%u_C ) ) deallocate( this%u_C  )
    
    if ( allocated(this%t_dn) ) deallocate( this%t_dn )
    if ( allocated(this%t_up) ) deallocate( this%t_up )
    
    if ( allocated(this%v_up) ) deallocate( this%v_up )
    if ( allocated(this%v_dn) ) deallocate( this%v_dn )
    
  end subroutine deallocate_boundaries_sub
  
end submodule Boundaries_init