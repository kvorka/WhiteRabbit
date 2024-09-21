submodule (Boundaries) Boundaries_initTemperature
  implicit none; contains
  
  module pure subroutine init_temp_up_sub(this)
    class(T_boundaries), intent(inout) :: this
    
    allocate( this%temp_up(this%jms) )
      this%temp_up = czero
    
  end subroutine init_temp_up_sub
  
end submodule Boundaries_initTemperature