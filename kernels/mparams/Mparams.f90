module Mparams
  use math
  
  type, public :: T_Mparams
    logical                        :: initvisc, visc_radial
    logical                        :: initlambda, lambda_radial
    logical                        :: initcp, cp_radial
    logical                        :: initalpha, alpha_radial
    integer                        :: nd, jmax, jms
    complex(kind=dbl), allocatable :: visc(:,:), lambda(:,:), cp(:,:), alpha(:,:)
    
    contains
    
    procedure :: init_sub       => init_mparams_sub
    procedure :: deallocate_sub => deallocate_mparams_sub

    procedure :: init_visc_sub
    procedure :: init_lambda_sub
    procedure :: init_cp_sub
    procedure :: init_alpha_sub
    
    procedure :: init_visc_radial_sub
    procedure :: init_lambda_radial_sub
    procedure :: init_cp_radial_sub
    procedure :: init_alpha_radial_sub
  
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
    
    module subroutine init_lambda_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_lambda_sub
    
    module subroutine init_cp_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_cp_sub
    
    module subroutine init_alpha_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_alpha_sub
    
    module subroutine init_visc_radial_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_visc_radial_sub
    
    module subroutine init_lambda_radial_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_lambda_radial_sub
    
    module subroutine init_cp_radial_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_cp_radial_sub
    
    module subroutine init_alpha_radial_sub(this)
      class(T_Mparams), intent(inout) :: this
    end subroutine init_alpha_radial_sub
  end interface

end module Mparams