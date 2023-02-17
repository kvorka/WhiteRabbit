module OceanMod
  use Math
  use PhysicalObject
  use OceanConstants
  implicit none

  type, extends(T_physicalObject), abstract, public :: T_ocean

    contains
    
    procedure :: init_ocean_sub       => init_ocean_sub
    procedure :: deallocate_ocean_sub => deallocate_ocean_sub

  end type T_ocean

  contains

  subroutine init_ocean_sub(this)
    class(T_ocean), intent(inout) :: this

    call this%init_objects_sub( nd = nd_ocean, jmax = jmax_ocean, r_ud = r_ud_ocean, &
                              & rgrid = grid_type_ocean, noharm = noharm_ocean       )
    
    call this%gravity%init_sub( gmod = gravity_ocean, g = 1 / this%ru**2 )

    this%n_iter = n_iter_ocean

    this%Pr = Pr_ocean
    this%Ra = Ra_ocean
    this%Ek = Ek_ocean

    this%D_ud         = D_ud_ocean
    this%rheology     = rheology_ocean
    this%mechanic_bnd = mechanic_bnd_ocean
    this%thermal_bnd  = thermal_bnd_ocean
    this%scaling      = scaling_ocean

  end subroutine init_ocean_sub

  subroutine deallocate_ocean_sub(this)
    class(T_ocean), intent(inout) :: this
    
    call this%deallocate_objects_sub()

  end subroutine deallocate_ocean_sub

end module OceanMod