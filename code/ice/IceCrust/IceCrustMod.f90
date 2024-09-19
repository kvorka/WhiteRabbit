module IceCrustMod
  use IceMod
  use IceTidesMod
  implicit none
  
  type, extends(T_ice), public :: T_iceCrust
    type(T_iceTides), private :: tides
    
    contains
    
    procedure :: init_sub        => init_iceCrust_sub
    procedure :: time_scheme_sub => iter_iceCrust_sub
    procedure :: solve_sub       => solve_iceCrust_sub

    procedure :: set_visc_sub    => visc_iceCrust_jm_sub
    procedure :: set_lambda_sub  => lambda_iceCrust_jm_sub
    
    procedure :: Vdelta_fn       => Vdelta_iceCrust_fn
    procedure :: set_layers_sub  => set_layers_iceCrust_sub
    
    procedure, private :: EE_sub => EE_iceCrust_sub
    
  end type T_iceCrust
  
  interface
    module subroutine init_iceCrust_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine init_iceCrust_sub
    
    module subroutine iter_iceCrust_sub(this, flux_bnd)
      class(T_iceCrust), intent(inout) :: this
      complex(kind=dbl), intent(in)    :: flux_bnd(:)
    end subroutine iter_iceCrust_sub
    
    module subroutine solve_iceCrust_sub(this, flux)
      class(T_iceCrust),           intent(inout) :: this
      complex(kind=dbl), optional, intent(in)    :: flux(:)
    end subroutine solve_iceCrust_sub
    
    module subroutine EE_iceCrust_sub(this, flux_bnd)
      class(T_iceCrust),           intent(inout) :: this
      complex(kind=dbl), optional, intent(in)    :: flux_bnd(:)
    end subroutine EE_iceCrust_sub
    
    module pure complex(kind=dbl) function Vdelta_iceCrust_fn(this, ir, ijm)
      class(T_iceCrust), intent(in) :: this
      integer,           intent(in) :: ir, ijm
    end function Vdelta_iceCrust_fn

    module subroutine set_layers_iceCrust_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine set_layers_iceCrust_sub
    
    module subroutine visc_iceCrust_jm_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine visc_iceCrust_jm_sub
    
    module subroutine lambda_iceCrust_jm_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine lambda_iceCrust_jm_sub
    
    module subroutine vypis_iceCrust_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine vypis_iceCrust_sub
  end interface
  
end module IceCrustMod