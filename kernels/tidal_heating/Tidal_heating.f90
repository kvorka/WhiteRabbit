module Tidal_heating
  use math
  implicit none
  
  type, public :: T_tidalHeating
    integer                        :: jms, nd
    complex(kind=dbl), allocatable :: htide(:,:)
    
    contains
    
    procedure, pass :: init_sub       => init_tidalHeating_sub
    procedure, pass :: deallocate_sub => deallocate_tidalHeating_sub
    
  end type T_tidalHeating
  
  interface
    module pure subroutine init_tidalHeating_sub(this, nd, jms)
      class(T_tidalHeating), intent(inout) :: this
      integer,               intent(in)    :: nd, jms
    end subroutine init_tidalHeating_sub
    
    module pure subroutine deallocate_tidalHeating_sub(this)
      class(T_tidalHeating), intent(inout) :: this
    end subroutine deallocate_tidalHeating_sub
  end interface
  
end module Tidal_heating