submodule (physicalobject) matrix_definitions
  implicit none ; contains
  
  module procedure mat_temp_fn
    
    select case (this%grid_type)
      case('homog')
        matica = matica_temp_hom_fn(this, j_in, a_in)
      
      case('chebv')
        matica = matica_temp_chb_fn(this, j_in, a_in)
    end select
    
  end procedure mat_temp_fn
  
  module procedure mat_torr_fn
    
    select case (this%scaling)
      case('christ')
        matica = matica_torr_chb_christ_viscos_fn(this, j_in, a_in)

      case default
        matica = matica_torr_chb_viscos_fn(this, j_in, a_in)
    end select
    
  end procedure mat_torr_fn
  
  module procedure mat_mech_fn
    
    select case(this%scaling)
      case('christ')
        matica = matica_mech_chb_christ_viscos_fn(this, j_in, a_in)

      case default 
        select case (this%grid_type)
          case('homog')
            select case (this%rheology)
              case('viscos')
                matica = matica_mech_hom_viscos_fn(this, j_in, a_in)
            
              case('viscel')
                matica = matica_mech_hom_viscel_fn(this, j_in, a_in)
            end select
      
          case('chebv')
            select case (this%rheology)
              case('viscos')
                matica = matica_mech_chb_viscos_fn(this, j_in, a_in)
            
              case('viscel')
                matica = matica_mech_chb_viscel_fn(this, j_in, a_in)
            end select
        end select
    end select
    
  end procedure mat_mech_fn
  
end submodule matrix_definitions