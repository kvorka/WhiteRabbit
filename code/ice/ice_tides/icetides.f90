module icetides
  use ice
  implicit none
    
  type, extends(T_ice), public :: T_iceTides
    
    contains
    
    procedure :: init_sub => init_iceTides_sub

    procedure :: compute_sub    => compute_iceTides_sub
    procedure :: Vdelta_fn      => Vdelta_iceTides_fn
    procedure :: set_layers_sub => set_layers_iceTides_sub
    
    procedure :: set_visc_sub => visc_iceTides_jm_sub
    
    procedure, private :: EE_mech_sub => EE_mech_iceTides_sub
    
  end type T_iceTides
  
  interface
    module subroutine init_iceTides_sub(this)
      class(T_iceTides), intent(inout) :: this
    end subroutine init_iceTides_sub
    
    module subroutine compute_iceTides_sub(this, visc_prof)
      class(T_iceTides), intent(inout) :: this
      complex(kind=dbl), intent(in)    :: visc_prof(:,:)
    end subroutine compute_iceTides_sub
    
    module pure complex(kind=dbl) function Vdelta_iceTides_fn(this, ir, ijm)
      class(T_iceTides), intent(in) :: this
      integer,           intent(in) :: ir, ijm
    end function Vdelta_iceTides_fn
    
    module subroutine set_layers_iceTides_sub(this)
      class(T_iceTides), intent(inout) :: this
    end subroutine set_layers_iceTides_sub
    
    module subroutine visc_iceTides_jm_sub(this)
      class(T_iceTides), intent(inout) :: this
    end subroutine visc_iceTides_jm_sub
    
    module subroutine EE_mech_iceTides_sub(this)
      class(T_iceTides), intent(inout) :: this
    end subroutine EE_mech_iceTides_sub
  end interface
  
end module icetides