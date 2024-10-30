submodule (mparams) cond
  implicit none; contains
  
  module procedure init_lambda_sub
    
    this%initlambda = .true.
    this%lambda_radial = .false.
    
    allocate( this%lambda(this%jms,this%nd) )
      this%lambda = cone
    
  end procedure init_lambda_sub
  
  module procedure init_lambda_radial_sub
    
    this%initlambda = .true.
    this%lambda_radial = .true.
    
    allocate( this%lambda(1,this%nd) )
      this%lambda = cone
    
  end procedure init_lambda_radial_sub
  
end submodule cond