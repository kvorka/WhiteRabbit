submodule (Mparams) Mparams_cap
  implicit none; contains
  
  module subroutine init_cap_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initcp = .true.
    
    allocate( this%cp(this%jms,this%nd) ) ; this%cp = czero
    
  end subroutine init_cap_sub
  
end submodule Mparams_cap