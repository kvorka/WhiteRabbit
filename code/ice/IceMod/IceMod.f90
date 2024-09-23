module iceMod
  use Lateral_integral
  use PhysicalObject
  use IceConstants
  use IceViscosity
  use IceConductivity
  use IceCapacity
  use IceExpansivity
  use IceBndTemp
  implicit none
  
  type, extends(T_physicalObject), abstract, public :: T_ice
    real(kind=dbl) :: diam, lambdaC, hC, lambdaU, viscU, cutoff, alphaU, cU, kappaU, Td, Tu, period, omega, g, mu
    real(kind=dbl) :: rC, rI2, rhoC, rhoI2, rhoW, rhoI
    
    contains
    
    procedure :: init_ice_sub      => init_ice_sub
    procedure :: deallocate_sub    => deallocate_ice_sub
    
    procedure :: avrg_temperature_ice_ir_fn, avrg_temperature_ice_irr_fn
    procedure :: avrg_stress_ice_ir_fn, avrg_stress_ice_irr_fn
    
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
    
    module pure real(kind=dbl) function avrg_temperature_ice_ir_fn(this, ir)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: ir
    end function avrg_temperature_ice_ir_fn
    
    module pure real(kind=dbl) function avrg_temperature_ice_irr_fn(this, ir)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: ir
    end function avrg_temperature_ice_irr_fn
    
    module pure real(kind=dbl) function avrg_stress_ice_ir_fn(this, ir)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: ir
    end function avrg_stress_ice_ir_fn
    
    module pure real(kind=dbl) function avrg_stress_ice_irr_fn(this, ir)
      class(T_ice),  intent(in) :: this
      integer,       intent(in) :: ir
    end function avrg_stress_ice_irr_fn
  end interface
  
end module iceMod
