submodule (boundaries) inittemp
  implicit none; contains
  
  module pure subroutine init_temp_up_sub(this)
    class(T_boundaries), intent(inout) :: this
    
    allocate( this%temp_up(this%jms) )
    
    call zero_carray_sub( this%jms, this%temp_up )
    
  end subroutine init_temp_up_sub
  
  module pure subroutine init_flux_up_sub(this)
    class(T_boundaries), intent(inout) :: this
    
    allocate( this%flux_up(this%jms) )
    
    call zero_carray_sub( this%jms, this%flux_up )
    
  end subroutine init_flux_up_sub
  
end submodule inittemp