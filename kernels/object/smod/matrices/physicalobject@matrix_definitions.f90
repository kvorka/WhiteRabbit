submodule (physicalobject) matrix_definitions
  implicit none ; contains
  
  module pure function mat_temp_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),        allocatable  :: matica(:,:)
    
    select case (this%grid_type)
      case('homog')
        matica = matica_temp_hom_fn(this, j_in, a_in)
      
      case('chebv')
        matica = matica_temp_chb_fn(this, j_in, a_in)
    end select
    
  end function mat_temp_fn
  
  module pure function mat_torr_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),        allocatable  :: matica(:,:)
    
    select case (this%scaling)
      case('christ')
        matica = matica_torr_chb_christ_viscos_fn(this, j_in, a_in)

      case default
        matica = matica_torr_chb_viscos_fn(this, j_in, a_in)
    end select
    
  end function mat_torr_fn
  
  module pure function mat_mech_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),        allocatable  :: matica(:,:)
    
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
    
  end function mat_mech_fn
  
end submodule matrix_definitions