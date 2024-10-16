module oceantides
  use ocean
  implicit none

  type, extends(T_ocean), public :: T_oceanTides
    complex(kind=dbl), allocatable, private :: v201(:), v203(:), v221(:), v223(:)
    real(kind=dbl)                          :: heating
    integer                                 :: number_of_periods, k_of_period

    contains

    procedure, public, pass :: init_sub        => init_oceanTides_sub
    procedure, public, pass :: iter_sub        => iter_oceanTides_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanTides_sub
    procedure, public, pass :: vypis_ocean_sub => vypis_oceanTides_sub

  end type T_oceanTides
  
  interface
    module subroutine init_oceanTides_sub(this)
      class(T_oceanTides), intent(inout) :: this
    end subroutine init_oceanTides_sub
    
    module subroutine time_scheme_oceanTides_sub(this)
      class(T_oceanTides), intent(inout) :: this
    end subroutine time_scheme_oceanTides_sub
    
    module subroutine iter_oceanTides_sub(this)
      class(T_oceanTides), intent(inout) :: this
    end subroutine iter_oceanTides_sub
    
    module subroutine vypis_oceanTides_sub(this)
      class(T_oceanTides), intent(inout) :: this
    end subroutine vypis_oceanTides_sub
  end interface
  
end module oceantides