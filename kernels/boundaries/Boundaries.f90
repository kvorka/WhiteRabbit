module Boundaries
  use math
  implicit none
  
  type, public :: T_boundaries
    integer                        :: jms
    complex(kind=dbl), allocatable :: temp_up(:), flux_up(:)
    complex(kind=dbl), allocatable :: u_dn(:), u_up(:), u_I2(:), u_C(:)
    complex(kind=dbl), allocatable :: t_dn(:), t_up(:)
    complex(kind=dbl), allocatable :: v_dn(:), v_up(:)
    
    contains
    
    procedure :: init_sub       => init_boundaries_sub
    procedure :: nulify_sub     => nulify_boundaries_sub
    procedure :: deallocate_sub => deallocate_boundaries_sub

    procedure :: init_layers_sub
    procedure :: init_layer_up_sub
    procedure :: init_temp_up_sub
    procedure :: init_flux_up_sub

  end type T_boundaries
  
  interface
    module pure subroutine init_boundaries_sub(this, jms)
      class(T_boundaries), intent(inout) :: this
      integer,             intent(in)    :: jms
    end subroutine init_boundaries_sub
    
    module pure subroutine deallocate_boundaries_sub(this)
      class(T_boundaries), intent(inout) :: this
    end subroutine deallocate_boundaries_sub
    
    module pure subroutine init_layers_sub(this)
      class(T_boundaries), intent(inout) :: this
    end subroutine init_layers_sub
    
    module pure subroutine init_layer_up_sub(this)
      class(T_boundaries), intent(inout) :: this
    end subroutine init_layer_up_sub
    
    module pure subroutine init_temp_up_sub(this)
      class(T_boundaries), intent(inout) :: this
    end subroutine init_temp_up_sub
    
    module pure subroutine init_flux_up_sub(this)
      class(T_boundaries), intent(inout) :: this
    end subroutine init_flux_up_sub
    
    module pure subroutine nulify_boundaries_sub(this)
      class(T_boundaries), intent(inout) :: this
    end subroutine nulify_boundaries_sub
  end interface
  
end module Boundaries