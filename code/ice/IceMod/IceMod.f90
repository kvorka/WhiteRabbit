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
    
    procedure :: htide_fn          => htide_ice_4_fn
    procedure :: tidal_heating_sub => tidal_heating_ice_4_sub
    
    procedure :: average_temperature_ice_ir_fn, average_temperature_ice_irr_fn
    procedure :: average_stress_ice_ir_fn, average_stress_ice_irr_fn
    
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
    
    module pure complex(kind=dbl) function htide_ice_4_fn(this, ir, ijm)
      class(T_ice), intent(in) :: this
      integer,      intent(in) :: ir, ijm
    end function htide_ice_4_fn

    module subroutine tidal_heating_ice_4_sub(this)
      class(T_ice), intent(inout) :: this
    end subroutine tidal_heating_ice_4_sub
    
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
