module iceMod
  use PhysicalObject
  use IceConstants
  use IceViscosity
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
    procedure :: lambda_fn         => lambda_ice_fn
    procedure :: cp_fn             => cp_ice_fn
    procedure :: alpha_fn          => alpha_ice_fn
    procedure :: visc_fn           => visc_ice_fn
    
    procedure :: temperature_ice_r_fn, temperature_ice_rr_fn, temperature_ice_r_jm_sub
    procedure :: devstress_ice_r_fn
    procedure :: visc_ice_jm_sub
    
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

    module pure real(kind=dbl) function lambda_ice_fn(this, i)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: i
    end function lambda_ice_fn
    
    module pure real(kind=dbl) function cp_ice_fn(this, i)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: i
    end function cp_ice_fn
    
    module pure real(kind=dbl) function alpha_ice_fn(this, i)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: i
    end function alpha_ice_fn
    
    module pure real(kind=dbl) function visc_ice_fn(this, i)
      class(T_ice),      intent(in) :: this
      integer,           intent(in) :: i
    end function visc_ice_fn
    
    module subroutine visc_ice_jm_sub(this)
      class(T_ice), intent(inout) :: this
    end subroutine visc_ice_jm_sub

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
    
    module pure real(kind=dbl) function temperature_ice_r_fn(this, i)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: i
    end function temperature_ice_r_fn
    
    module pure subroutine temperature_ice_r_jm_sub(this, i, temperature)
      class(T_ice),      intent(in)  :: this
      integer,           intent(in)  :: i
      complex(kind=dbl), intent(out) :: temperature(:)
    end subroutine temperature_ice_r_jm_sub
    
    module pure real(kind=dbl) function temperature_ice_rr_fn(this, i)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: i
    end function temperature_ice_rr_fn
    
    module pure real(kind=dbl) function devstress_ice_r_fn(this, i)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: i
    end function devstress_ice_r_fn
  end interface
  
end module iceMod
