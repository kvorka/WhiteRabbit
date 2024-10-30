submodule (physicalobject) coriolis
  implicit none; contains
  
  module procedure coriolis_rr_jml_sub
    
    select case (this%scaling)
      case ('christ')
        call ezvv_sub(this%jmax, 2._dbl, v, coriolis)
      
      case('basics')
        call ezvv_sub(this%jmax, 2/this%Ek, v, coriolis)
    end select
    
  end procedure coriolis_rr_jml_sub
  
end submodule coriolis