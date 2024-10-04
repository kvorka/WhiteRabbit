module oceanconv
  use ocean
  implicit none
  
  type, extends(T_ocean), public :: T_oceanConv
    contains
    
    procedure, public, pass :: init_sub        => init_oceanConv_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanConv_sub
  end type T_oceanConv
  
  interface
    module subroutine init_oceanConv_sub(this)
      class(T_oceanConv), intent(inout) :: this
    end subroutine init_oceanConv_sub
    
    module subroutine time_scheme_oceanConv_sub(this)
      class(T_oceanConv), intent(inout) :: this
    end subroutine time_scheme_oceanConv_sub
  end interface
  
end module oceanconv
