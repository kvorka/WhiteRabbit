module Mparams
  use Math
  
  type, public :: T_Mparams
    logical                        :: initvisc, initviscradial
    logical                        :: initlambda, initcp
    integer                        :: nd, jmax, jms
    complex(kind=dbl), allocatable :: visc_radial(:)
    complex(kind=dbl), allocatable :: visc(:,:), lambda(:,:), cp(:,:)
    
    contains
    
    procedure :: init_sub              => init_mparams_sub
    procedure :: deallocate_sub        => deallocate_mparams_sub

    procedure :: init_visc_sub
    procedure :: init_conductivity_sub
    procedure :: init_capacity_sub
    
    procedure :: init_visc_radial_sub
  
  end type T_Mparams
  
  interface
    module subroutine init_mparams_sub(this, nd, jmax)
      class(T_Mparams), intent(inout) :: this
      integer,          intent(in)    :: nd, jmax
    end subroutine init_mparams_sub
    
    module subroutine deallocate_mparams_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine deallocate_mparams_sub
    
    module subroutine init_visc_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_visc_sub
    
    module subroutine init_visc_radial_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_visc_radial_sub
    
    module subroutine init_cond_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_cond_sub
    
    module subroutine init_cap_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_cap_sub
  end interface

end module Mparams