submodule (Mparams) Mparams_cap
  implicit none; contains
  
  module subroutine init_cp_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initcp = .true.
    this%cp_radial = .false.
    
    allocate( this%cp(this%jms,this%nd+1) )
      this%cp = cone
    
  end subroutine init_cp_sub
  
  module subroutine init_cp_radial_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initcp = .true.
    this%cp_radial = .true.
    
    allocate( this%cp(1,this%nd+1) )
      this%cp = cone
    
  end subroutine init_cp_radial_sub
  
end submodule Mparams_cap