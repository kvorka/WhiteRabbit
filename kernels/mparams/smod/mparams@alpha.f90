submodule (mparams) alpha
  implicit none; contains
  
  module subroutine init_alpha_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initalpha = .true.
    this%alpha_radial = .false.
    
    allocate( this%alpha(this%jms,this%nd+1) )
      this%alpha = cone
    
  end subroutine init_alpha_sub
  
  module subroutine init_alpha_radial_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initalpha = .true.
    this%alpha_radial = .true.
    
    allocate( this%alpha(1,this%nd+1) )
      this%alpha = cone
    
  end subroutine init_alpha_radial_sub
  
end submodule alpha