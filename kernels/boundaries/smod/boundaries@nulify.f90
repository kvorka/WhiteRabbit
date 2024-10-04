submodule (boundaries) nulify
  implicit none; contains
  
  module pure subroutine nulify_boundaries_sub(this)
    class(T_boundaries), intent(inout) :: this
    
    if ( allocated(this%temp_up) ) call zero_carray_sub( this%jms, this%temp_up)
    if ( allocated(this%flux_up) ) call zero_carray_sub( this%jms, this%flux_up)
    
    if ( allocated(this%u_dn) ) call zero_carray_sub( this%jms, this%u_dn )
    if ( allocated(this%u_up) ) call zero_carray_sub( this%jms, this%u_up )
    if ( allocated(this%u_I2) ) call zero_carray_sub( this%jms, this%u_I2 )
    if ( allocated(this%u_C ) ) call zero_carray_sub( this%jms, this%u_C  )
    
    if ( allocated(this%t_dn) ) call zero_carray_sub( this%jms, this%t_dn )
    if ( allocated(this%t_up) ) call zero_carray_sub( this%jms, this%t_up )
    
    if ( allocated(this%v_up) ) call zero_carray_sub( this%jms, this%v_up )
    if ( allocated(this%v_dn) ) call zero_carray_sub( this%jms, this%v_dn )
    
  end subroutine nulify_boundaries_sub
  
end submodule nulify