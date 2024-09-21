submodule (Boundaries) Boundaries_nulify
  implicit none; contains
  
  module pure subroutine nulify_boundaries_sub(this)
    class(T_boundaries), intent(inout) :: this
    
    if ( allocated(this%temp_up) ) this%temp_up = czero
    
    if ( allocated(this%u_dn) ) this%u_dn = czero
    if ( allocated(this%u_up) ) this%u_up = czero
    if ( allocated(this%u_I2) ) this%u_I2 = czero
    if ( allocated(this%u_C ) ) this%u_C  = czero
    
    if ( allocated(this%t_dn) ) this%t_dn = czero
    if ( allocated(this%t_up) ) this%t_up = czero
    
    if ( allocated(this%v_up) ) this%v_up = czero
    if ( allocated(this%v_dn) ) this%v_dn = czero
    
  end subroutine nulify_boundaries_sub
  
end submodule Boundaries_nulify