submodule (oceantides) output
  implicit none; contains
  
  module procedure vypis_oceanTides_sub

    write(11,*) this%number_of_periods, stress_dim * this%heating / 1e6
    
  end procedure vypis_oceanTides_sub
  
end submodule output