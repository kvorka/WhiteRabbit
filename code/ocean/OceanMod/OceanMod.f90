module OceanMod
  use physicalobject
  use OceanConstants
  implicit none

  type, extends(T_physicalObject), abstract, public :: T_ocean
    complex(kind=dbl), allocatable :: nsph1(:,:), nsph2(:,:), ntorr(:,:), ntemp(:,:)
    
    contains
    
    procedure :: init_ocean_sub
    procedure :: speed_sub
    procedure :: set_boundary_deformation_sub
    
    procedure :: vypis_ocean_sub => vypis_ocean_sub
    procedure :: iter_sub        => iter_ocean_sub
    procedure :: init_state_sub  => init_state_ocean_sub

    procedure :: fullnl_sub          => fullnl_ocean_sub
    procedure :: coriolis_sub        => coriolis_ocean_sub
    procedure :: coriolis_vgradv_sub => coriolis_vgradv_ocean_sub

    procedure :: deallocate_sub  => deallocate_ocean_sub

    procedure(time_scheme_abstract), deferred, pass :: time_scheme_sub

  end type T_ocean

  abstract interface
    subroutine time_scheme_abstract(this)
       import :: T_ocean, dbl
       class(T_ocean), intent(inout) :: this
    end subroutine time_scheme_abstract
  end interface
  
  interface
    module subroutine init_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine init_ocean_sub
    
    module subroutine deallocate_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine deallocate_ocean_sub
    
    module subroutine init_state_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine init_state_ocean_sub
    
    module subroutine iter_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine iter_ocean_sub
    
    module subroutine speed_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine speed_sub
    
    module subroutine vypis_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine vypis_ocean_sub
    
    module subroutine set_boundary_deformation_sub(this, u_up, t_up)
      class(T_ocean),    intent(inout) :: this
      complex(kind=dbl), intent(in)    :: u_up(:), t_up(:)
    end subroutine set_boundary_deformation_sub
    
    module pure subroutine coriolis_ocean_sub(this, i)
      class(T_ocean), intent(inout) :: this
      integer,        intent(in)    :: i
    end subroutine coriolis_ocean_sub
    
    module subroutine coriolis_vgradv_ocean_sub(this, i)
      class(T_ocean),    intent(inout) :: this
      integer,           intent(in)    :: i
    end subroutine coriolis_vgradv_ocean_sub
    
    module subroutine fullnl_ocean_sub(this, i)
      class(T_ocean), intent(inout) :: this
      integer,        intent(in)    :: i
    end subroutine fullnl_ocean_sub
  end interface
  
end module OceanMod