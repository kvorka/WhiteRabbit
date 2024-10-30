submodule (mparams) visc
  implicit none; contains
  
  module procedure init_visc_sub
    
    this%initvisc = .true.
    this%visc_radial = .false.
    
    allocate( this%visc(this%jms,this%nd) )
      this%visc = cone
    
  end procedure init_visc_sub
  
  module procedure init_visc_radial_sub
    
    this%initvisc = .true.
    this%visc_radial = .true.
    
    allocate( this%visc(1,this%nd) )
      this%visc = cone
    
  end procedure init_visc_radial_sub
  
end submodule visc