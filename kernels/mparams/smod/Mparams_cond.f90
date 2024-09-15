submodule (Mparams) Mparams_cond
  implicit none; contains
  
  module subroutine init_cond_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initlambda = .true.
    
    allocate( this%lambda(this%jms,this%nd) ) ; this%lambda = czero
    
  end subroutine init_cond_sub
  
end submodule Mparams_cond