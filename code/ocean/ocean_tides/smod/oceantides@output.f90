submodule (oceantides) output
  implicit none; contains
  
  module subroutine vypis_oceanTides_sub(this)
    class(T_oceanTides), intent(inout) :: this

    write(11,*) this%number_of_periods, stress_dim * this%heating / 1e6
    
  end subroutine vypis_oceanTides_sub
  
end submodule output