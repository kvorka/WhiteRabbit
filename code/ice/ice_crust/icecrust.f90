module icecrust
  use ice
  use icetides
  implicit none
  
  type, extends(T_ice), public :: T_iceCrust
    type(T_iceTides), private :: tides
    
    contains
    
    procedure :: init_sub        => init_iceCrust_sub
    procedure :: time_scheme_sub => iter_iceCrust_sub

    procedure :: solve_sub            => solve_iceCrust_sub
    procedure :: solve_conduction_sub => solve_conduction_iceCrust_sub
    
    procedure :: set_surfTemp_sub => surfTemp_iceCrust_jm_sub
    procedure :: set_lambda_sub   => lambda_iceCrust_jm_sub
    procedure :: set_cp_sub       => cp_iceCrust_jm_sub
    procedure :: set_alpha_sub    => alpha_iceCrust_jm_sub
    procedure :: set_visc_sub     => visc_iceCrust_jm_sub
    
    procedure :: Vdelta_fn       => Vdelta_iceCrust_fn
    procedure :: set_layers_sub  => set_layers_iceCrust_sub
    
    procedure :: vypis_iceCrust_sub
    
    procedure :: cpdivq_sub
    procedure :: mvgradT_sub
    procedure :: mvgradT_cpdivq_sub
    procedure :: mlambdagradT_sub
    
    procedure :: EE_sub      => EE_iceCrust_sub
    procedure :: EE_temp_sub => EE_temp_iceCrust_sub
    procedure :: EE_mech_sub => EE_mech_iceCrust_sub
    
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
      class(T_iceCrust), intent(inout) :: this
      complex(kind=dbl), intent(in)    :: flux(:)
    end subroutine solve_iceCrust_sub
    
    module subroutine solve_conduction_iceCrust_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine solve_conduction_iceCrust_sub
    
    module subroutine EE_iceCrust_sub(this, flux_bnd)
      class(T_iceCrust),           intent(inout) :: this
      complex(kind=dbl), optional, intent(in)    :: flux_bnd(:)
    end subroutine EE_iceCrust_sub
    
    module subroutine EE_temp_iceCrust_sub(this, flux)
      class(T_iceCrust), intent(inout) :: this
      complex(kind=dbl), intent(inout) :: flux(:)
    end subroutine EE_temp_iceCrust_sub
    
    module subroutine EE_mech_iceCrust_sub(this, flux)
      class(T_iceCrust), intent(inout) :: this
      complex(kind=dbl), intent(in)    :: flux(:)
    end subroutine EE_mech_iceCrust_sub
    
    module pure complex(kind=dbl) function Vdelta_iceCrust_fn(this, ir, ijm)
      class(T_iceCrust), intent(in) :: this
      integer,           intent(in) :: ir, ijm
    end function Vdelta_iceCrust_fn

    module subroutine set_layers_iceCrust_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine set_layers_iceCrust_sub
    
    module subroutine mvgradT_sub(this, ir, mvgradT)
      class(T_iceCrust), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: mvgradT(:)
    end subroutine mvgradT_sub
    
    module subroutine mvgradT_cpdivq_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine mvgradT_cpdivq_sub
    
    module subroutine mlambdagradT_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine mlambdagradT_sub
    
    module subroutine cpdivq_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine cpdivq_sub
    
    module subroutine surfTemp_iceCrust_jm_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine surfTemp_iceCrust_jm_sub
    
    module subroutine visc_iceCrust_jm_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine visc_iceCrust_jm_sub
    
    module subroutine lambda_iceCrust_jm_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine lambda_iceCrust_jm_sub
    
    module subroutine cp_iceCrust_jm_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine cp_iceCrust_jm_sub
    
    module subroutine alpha_iceCrust_jm_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine alpha_iceCrust_jm_sub
    
    module subroutine vypis_iceCrust_sub(this)
      class(T_iceCrust), intent(inout) :: this
    end subroutine vypis_iceCrust_sub
  end interface
  
end module icecrust