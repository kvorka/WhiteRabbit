module MatrixDefinitions
  use ThermalMatrices
  use ToroidalVisc
  use SpheroidalVisc
  use SpheroidalViscel
  implicit none
  
  public :: matica_temp_fn
  public :: matica_torr_fn
  public :: matica_mech_fn
  
  contains
  
  function matica_temp_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),        allocatable  :: matica(:,:)
    
    select case (this%grid_type)
      case('homog')
        allocate( matica(7,3*this%nd+1) ); matica = matica_temp_hom_fn(this, j_in, a_in)
      
      case('chebv')
        allocate(matica(11,3*this%nd+1) ); matica = matica_temp_chb_fn(this, j_in, a_in)
    end select
    
  end function matica_temp_fn
  
  function matica_torr_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),        allocatable  :: matica(:,:)
    
    select case (this%scaling)
      case('christ')
        allocate( matica(11, 3*this%nd+1) ); matica = matica_torr_chb_christ_viscos_fn(this, j_in, a_in)

      case default
        allocate( matica(11, 3*this%nd+1) ); matica = matica_torr_chb_viscos_fn(this, j_in, a_in)

    end select
    
  end function matica_torr_fn
  
  function matica_mech_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),        allocatable  :: matica(:,:)
    
    select case(this%scaling)
      case('christ')
        allocate( matica(23,6*this%nd+2) ); matica = matica_mech_chb_christ_viscos_fn(this, j_in, a_in)

      case default 
        select case (this%grid_type)
          case('homog')
            allocate( matica(15,6*this%nd+2) )
        
            select case (this%rheology)
              case('viscos')
                matica = matica_mech_hom_viscos_fn(this, j_in, a_in)
            
              case('viscel')
                matica = matica_mech_hom_viscel_fn(this, j_in, a_in)
            end select
      
          case('chebv')
            allocate( matica(23,6*this%nd+2) )
        
            select case (this%rheology)
              case('viscos')
                matica = matica_mech_chb_viscos_fn(this, j_in, a_in)
            
              case('viscel')
                matica = matica_mech_chb_viscel_fn(this, j_in, a_in)
            end select
        end select
    end select
    
  end function matica_mech_fn
  
end module MatrixDefinitions