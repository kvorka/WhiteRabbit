submodule (mparams) visc
  implicit none; contains
  
  module subroutine init_visc_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initvisc = .true.
    this%visc_radial = .false.
    
    allocate( this%visc(this%jms,this%nd) )
      this%visc = cone
    
  end subroutine init_visc_sub
  
  module subroutine init_visc_radial_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initvisc = .true.
    this%visc_radial = .true.
    
    allocate( this%visc(1,this%nd) )
      this%visc = cone
    
  end subroutine init_visc_radial_sub

end submodule visc