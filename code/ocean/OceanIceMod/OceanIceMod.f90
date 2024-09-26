module OceanIceMod
  use OceanMod
  implicit none
  
  type, extends(T_ocean), public :: T_oceanice
    contains
    
    procedure, public, pass :: init_sub        => init_oceanice_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanice_sub
  end type T_oceanice
  
  interface
    module subroutine init_oceanice_sub(this)
      class(T_oceanice), intent(inout) :: this
    end subroutine init_oceanice_sub
    
    module subroutine time_scheme_oceanice_sub(this)
      class(T_oceanIce), intent(inout) :: this
    end subroutine time_scheme_oceanice_sub
  end interface
  
end module OceanIceMod