module oceanice
  use ocean
  implicit none
  
  type, extends(T_ocean), public :: T_oceanice
    contains
    
    procedure, public, pass :: init_sub        => init_oceanice_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanice_sub
    procedure, public, pass :: iter_sub        => iter_oceanIce_sub
    procedure, public, pass :: set_ubnd_sub    => set_ubnd_oceanIce_sub
  end type T_oceanice
  
  interface
    module subroutine init_oceanice_sub(this)
      class(T_oceanice), intent(inout) :: this
    end subroutine init_oceanice_sub
    
    module subroutine set_ubnd_oceanIce_sub(this, u_up, t_up)
      class(T_oceanice), intent(inout) :: this
      complex(kind=dbl), intent(in)    :: u_up(:), t_up(:)
    end subroutine set_ubnd_oceanIce_sub
    
    module subroutine time_scheme_oceanice_sub(this)
      class(T_oceanIce), intent(inout) :: this
    end subroutine time_scheme_oceanice_sub
    
    module subroutine iter_oceanIce_sub(this)
      class(T_oceanice), intent(inout) :: this
    end subroutine iter_oceanIce_sub
  end interface
  
end module oceanice