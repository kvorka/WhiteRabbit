submodule (mparams) cap
  implicit none; contains
  
  module procedure init_cp_sub
    
    this%initcp = .true.
    this%cp_radial = .false.
    
    allocate( this%cp(this%jms,this%nd+1) )
      this%cp = cone
    
  end procedure init_cp_sub
  
  module procedure init_cp_radial_sub
    
    this%initcp = .true.
    this%cp_radial = .true.
    
    allocate( this%cp(1,this%nd+1) )
      this%cp = cone
    
  end procedure init_cp_radial_sub
  
end submodule cap