module iceMod
  use PhysicalObject
  use IceConstants
  use IceViscosity
  use IceConductivity
  use IceCapacity
  use IceExpansivity
  implicit none
  
  type, extends(T_physicalObject), abstract, public :: T_ice
    logical        :: andrade
    real(kind=dbl) :: diam, lambdaC, hC, lambdaU, viscU, cutoff, alphaU, cU, kappaU, Td, Tu, period, omega, g, mu
    real(kind=dbl) :: rC, rI2, rhoC, rhoI2, rhoW, rhoI
    
    contains
    
    procedure :: init_ice_sub      => init_ice_sub
    procedure :: deallocate_sub    => deallocate_ice_sub
    
    procedure :: Vdelta_fn         => Vdelta_ice_fn
    procedure :: htide_fn          => htide_ice_4_fn
    procedure :: set_layers_sub    => set_layers_ice_sub
    procedure :: tidal_heating_sub => tidal_heating_ice_4_sub
    
    procedure :: lambda_r_fn  => lambda_r_ice_fn
    procedure :: lambda_rr_fn => lambda_rr_ice_fn
    procedure :: cp_r_fn      => cp_r_ice_fn
    procedure :: cp_rr_fn     => cp_rr_ice_fn
    procedure :: alpha_r_fn   => alpha_r_ice_fn
    procedure :: alpha_rr_fn  => alpha_rr_ice_fn
    procedure :: visc_r_fn    => visc_r_ice_fn
    procedure :: visc_rr_fn   => visc_rr_ice_fn
    
    procedure :: average_temperature_ice_ir_fn, average_temperature_ice_irr_fn
    procedure :: average_stress_ice_ir_fn, average_stress_ice_irr_fn
    procedure :: visc_ice_jm_sub, lambda_ice_jm_sub
    
  end type T_ice

  interface
    module subroutine init_ice_sub(this, jmax_in, rheol_in, n_iter, noharm)
      class(T_ice),      intent(inout) :: this
      integer,           intent(in)    :: jmax_in, n_iter
      character(len=*),  intent(in)    :: rheol_in
      logical, optional, intent(in)    :: noharm
    end subroutine init_ice_sub
    
    module subroutine deallocate_ice_sub(this)
      class(T_ice), intent(inout) :: this
    end subroutine deallocate_ice_sub

    module pure real(kind=dbl) function lambda_r_ice_fn(this, ir)
      class(T_ice), intent(in) :: this
      integer,      intent(in) :: ir
    end function lambda_r_ice_fn
    
    module pure real(kind=dbl) function lambda_rr_ice_fn(this, ir)
      class(T_ice), intent(in) :: this
      integer,      intent(in) :: ir
    end function lambda_rr_ice_fn
    
    module pure real(kind=dbl) function cp_r_ice_fn(this, ir)
      class(T_ice), intent(in) :: this
      integer,      intent(in) :: ir
    end function cp_r_ice_fn
    
    module pure real(kind=dbl) function cp_rr_ice_fn(this, ir)
      class(T_ice), intent(in) :: this
      integer,      intent(in) :: ir
    end function cp_rr_ice_fn
    
    module pure real(kind=dbl) function alpha_r_ice_fn(this, ir)
      class(T_ice), intent(in) :: this
      integer,      intent(in) :: ir
    end function alpha_r_ice_fn
    
    module pure real(kind=dbl) function alpha_rr_ice_fn(this, ir)
      class(T_ice), intent(in) :: this
      integer,      intent(in) :: ir
    end function alpha_rr_ice_fn
    
    module pure real(kind=dbl) function visc_r_ice_fn(this, ir)
      class(T_ice), intent(in) :: this
      integer,      intent(in) :: ir
    end function visc_r_ice_fn
    
    module pure real(kind=dbl) function visc_rr_ice_fn(this, ir)
      class(T_ice), intent(in) :: this
      integer,      intent(in) :: ir
    end function visc_rr_ice_fn
    
    module subroutine visc_ice_jm_sub(this)
      class(T_ice), intent(inout) :: this
    end subroutine visc_ice_jm_sub
    
    module subroutine lambda_ice_jm_sub(this)
      class(T_ice), intent(inout) :: this
    end subroutine lambda_ice_jm_sub

    module pure complex(kind=dbl) function htide_ice_4_fn(this, ir, ijm)
      class(T_ice), intent(in) :: this
      integer,      intent(in) :: ir, ijm
    end function htide_ice_4_fn

    module subroutine tidal_heating_ice_4_sub(this)
      class(T_ice), intent(inout) :: this
    end subroutine tidal_heating_ice_4_sub

    module pure complex(kind=dbl) function Vdelta_ice_fn(this, ir, ijm)
      class(T_ice),       intent(in) :: this
      integer,            intent(in) :: ir, ijm
    end function Vdelta_ice_fn

    module subroutine set_layers_ice_sub(this)
      class(T_ice),      intent(inout) :: this
    end subroutine set_layers_ice_sub
    
    module pure real(kind=dbl) function average_temperature_ice_ir_fn(this, ir)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: ir
    end function average_temperature_ice_ir_fn
    
    module pure real(kind=dbl) function average_temperature_ice_irr_fn(this, ir)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: ir
    end function average_temperature_ice_irr_fn
    
    module pure real(kind=dbl) function average_stress_ice_ir_fn(this, ir)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: ir
    end function average_stress_ice_ir_fn
    
    module pure real(kind=dbl) function average_stress_ice_irr_fn(this, ir)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: ir
    end function average_stress_ice_irr_fn
  end interface
  
end module iceMod
