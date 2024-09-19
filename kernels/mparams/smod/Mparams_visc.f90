submodule (Mparams) Mparams_visc
  implicit none; contains
  
  module subroutine init_visc_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initvisc = .true.
    
    allocate( this%visc(this%jms,this%nd) ) ; this%visc = czero
    
  end subroutine init_visc_sub
  
  module subroutine init_visc_radial_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initviscradial = .true.
    
    allocate( this%visc(this%nd) ) ; this%visc_radial = zero
    
  end subroutine init_visc_radial_sub

end submodule Mparams_visc