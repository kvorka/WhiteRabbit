submodule (Mparams) Mparams_visc
  implicit none; contains
  
  module subroutine init_visc_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    this%initvisc = .true.
    
    allocate( this%visc(this%jms,this%nd) ) ; this%visc = czero
    
  end subroutine init_visc_sub

end submodule Mparams_visc