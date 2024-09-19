submodule (Mparams) Mparams_cond
  implicit none; contains
  
  module subroutine init_lambda_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initlambda = .true.
    this%lambda_radial = .false.
    
    allocate( this%lambda(this%jms,this%nd) )
      this%lambda = cone
    
  end subroutine init_lambda_sub
  
  module subroutine init_lambda_radial_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initlambda = .true.
    this%lambda_radial = .true.
    
    allocate( this%lambda(1,this%nd) )
      this%lambda = cone
    
  end subroutine init_lambda_radial_sub
  
end submodule Mparams_cond