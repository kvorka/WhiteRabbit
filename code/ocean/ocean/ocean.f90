module ocean
  use physicalobject
  use ocean_constants
  use omp_lib
  implicit none

  type, extends(T_physicalObject), abstract, public :: T_ocean
    complex(kind=dbl), allocatable :: nsph1(:,:), nsph2(:,:), ntorr(:,:), ntemp(:,:)
    
    contains
    
    procedure, public, pass :: init_ocean_sub
    
    
    procedure, public, pass :: vypis_ocean_sub    => vypis_ocean_sub
    procedure, public, pass :: init_state_sub     => init_state_ocean_sub
    procedure, public, pass :: init_temp_bbnd_sub => init_temp_bbnd_ocean_sub
    procedure, public, pass :: fullnl_sub         => fullnl_ocean_sub
    
    procedure, public, pass :: deallocate_sub  => deallocate_ocean_sub
    
    procedure(time_scheme_abstract), deferred, public, pass :: time_scheme_sub
    procedure(iter_sub),             deferred, public, pass :: iter_sub
    
  end type T_ocean
  
  abstract interface
    subroutine time_scheme_abstract(this)
       import :: T_ocean
       class(T_ocean), intent(inout) :: this
    end subroutine time_scheme_abstract
    
    subroutine iter_sub(this)
      import :: T_ocean
      class(T_ocean), intent(inout) :: this
    end subroutine iter_sub
  end interface
  
  interface
    module subroutine init_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine init_ocean_sub
    
    module subroutine deallocate_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine deallocate_ocean_sub
    
    module subroutine init_temp_bbnd_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine init_temp_bbnd_ocean_sub
    
    module subroutine init_state_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine init_state_ocean_sub
    
    module subroutine vypis_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine vypis_ocean_sub
    
    module subroutine fullnl_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine fullnl_ocean_sub
  end interface
  
end module ocean