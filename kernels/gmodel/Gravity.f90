module Gravity
  use conversions
  use RadialGrid
  implicit none
  
  real(kind=dbl), parameter :: kappa = 6.670d-11
  
  type, public :: T_gravity
    character(len=3), private :: gmod
    real(kind=dbl),   private :: Dcrust, g, exc, omega
    
    contains
    
    procedure :: init_sub       => init_gravity_sub
    procedure :: set_sub        => set_gravity_sub
    procedure :: deallocate_sub => deallocate_gravity_sub
    procedure :: g_fn, V_tide_fn, V_bnd_fn, V_rho_fn, V_rt_fn
    
  end type T_gravity
  
  interface
    module pure subroutine init_gravity_sub(this, gmod, g)
      class(T_gravity), intent(inout) :: this
      character(len=*), intent(in)    :: gmod
      real(kind=dbl),   intent(in)    :: g
    end subroutine init_gravity_sub
    
    module subroutine set_gravity_sub(this, Dcrust, omega, exc)
      class(T_gravity),           intent(inout) :: this
      real(kind=dbl),   optional, intent(in)    :: Dcrust, omega, exc
    end subroutine set_gravity_sub

    module pure real(kind=dbl) function g_fn(this, ri)
      class(T_gravity), intent(in) :: this
      real(kind=dbl),   intent(in) :: ri
    end function g_fn
    
    module pure complex(kind=dbl) function V_bnd_fn(this, j, m, ri, rb, drho, ujm)
      class(T_gravity),  intent(in) :: this
      integer,           intent(in) :: j, m
      real(kind=dbl),    intent(in) :: ri, rb, drho
      complex(kind=dbl), intent(in) :: ujm
    end function V_bnd_fn
    
    module pure complex(kind=dbl) function V_rho_fn(this, j, m, ri, field, rad_grid)
      class(T_gravity),    intent(in) :: this
      integer,             intent(in) :: j, m
      real(kind=dbl),      intent(in) :: ri
      complex(kind=dbl),   intent(in) :: field(:)
      class(T_radialGrid), intent(in) :: rad_grid
    end function V_rho_fn
    
    module pure complex(kind=dbl) function V_tide_fn(this, j, m, ri, phase)
      class(T_gravity), intent(in) :: this
      integer,          intent(in) :: j, m
      real(kind=dbl),   intent(in) :: ri, phase
    end function V_tide_fn
    
    module pure complex(kind=dbl) function V_rt_fn(this, j, m, ri)
      class(T_gravity), intent(in) :: this
      integer,          intent(in) :: j, m
      real(kind=dbl),   intent(in) :: ri
    end function V_rt_fn
    
    module pure subroutine deallocate_gravity_sub(this)
      class(T_gravity), intent(inout) :: this
    end subroutine deallocate_gravity_sub
  end interface
  
end module Gravity