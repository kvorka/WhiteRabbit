module Solution
  use Math
  implicit none

  type, public :: T_solution 
    integer                        :: nd, jmax, jms, jmv, jmt
    complex(kind=dbl), allocatable :: u_dn(:), u_up(:), u_I2(:), u_C(:), v_dn(:), v_up(:), t_dn(:), t_up(:), &
                                    & mech(:,:), temp(:,:), torr(:,:)
    
    contains
    
    procedure :: init_sub       => init_solution_sub
    procedure :: deallocate_sub => deallocate_solution_sub

    procedure :: init_stemp_sub
    procedure :: init_storr_sub
    procedure :: init_smech_sub
    procedure :: init_layers_sub
    procedure :: init_layer_u_sub
    procedure :: temp_fn
    procedure :: flux_fn
    procedure :: velocity_fn
    procedure :: deviatoric_stress_fn
    procedure :: temp_i_fn
    procedure :: flux_i_fn
    procedure :: velocity_i_fn
    procedure :: deviatoric_stress_i_fn
    procedure :: temp_jm_fn
    procedure :: flux_jml_fn
    procedure :: velocity_jml_fn
    procedure :: conv_velocity_jml_fn
    procedure :: deviatoric_stress_jml2_fn
    
  end type T_solution

  private :: init_solution_sub
  private :: deallocate_solution_sub

  interface
    module subroutine init_stemp_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_stemp_sub

    module subroutine init_storr_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_storr_sub

    module subroutine init_smech_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_smech_sub

    module subroutine init_layers_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_layers_sub

    module subroutine init_layer_u_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_layer_u_sub
  end interface

  interface
    module pure complex(kind=dbl) function temp_fn(this, i, j, m)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: i, j, m
    end function temp_fn

    module pure complex(kind=dbl) function flux_fn(this, i, j, m, l)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: i, j, m, l
    end function flux_fn

    module pure complex(kind=dbl) function velocity_fn(this, i, j, m, l)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: i, j, m, l
    end function velocity_fn

    module pure complex(kind=dbl) function deviatoric_stress_fn(this, i, j, m, l)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: i, j, m, l
    end function deviatoric_stress_fn
  end interface

  interface
    module pure function temp_i_fn(this, j, m) result(temp)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: j, m
      complex(kind=dbl)             :: temp(this%nd+1)
    end function temp_i_fn

    module pure function flux_i_fn(this, j, m, l) result(flux)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: j, m, l
      complex(kind=dbl)             :: flux(this%nd)
    end function flux_i_fn

    module pure function velocity_i_fn(this, j, m, l) result(velocity)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: j, m, l
      complex(kind=dbl)             :: velocity(this%nd+1)
    end function velocity_i_fn

    module pure function deviatoric_stress_i_fn(this, j, m, l) result(dstress)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: j, m, l
      complex(kind=dbl)             :: dstress(this%nd)
    end function deviatoric_stress_i_fn
  end interface

  interface
    module pure function temp_jm_fn(this, i) result(temp)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: i
      complex(kind=dbl)             :: temp(this%jms)
    end function temp_jm_fn

    module pure function flux_jml_fn(this, i) result(flux)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: i
      complex(kind=dbl)             :: flux(this%jmv)
    end function flux_jml_fn

    module pure function velocity_jml_fn(this, i) result(velocity)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: i
      complex(kind=dbl)             :: velocity(this%jmv)
    end function velocity_jml_fn

    module pure function conv_velocity_jml_fn(this, i) result(velocity)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: i
      complex(kind=dbl)             :: velocity(this%jmv)
    end function conv_velocity_jml_fn

    module pure function deviatoric_stress_jml2_fn(this, i) result(dstress)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: i
      complex(kind=dbl)             :: dstress(this%jmt)
    end function deviatoric_stress_jml2_fn
  end interface

  contains

  subroutine init_solution_sub(this, nd, jmax)
    class(T_solution), intent(inout) :: this
    integer,           intent(in)    :: nd, jmax
    
    this%nd   = nd
    this%jmax = jmax
    this%jms  =     jmax * (jmax+1) / 2 + jmax   + 1
    this%jmv  = 3*( jmax * (jmax+1) / 2 + jmax ) + 1
    this%jmt  = 5*( jmax * (jmax+1) / 2 + jmax ) + 1

  end subroutine init_solution_sub

  subroutine deallocate_solution_sub(this)
    class(T_solution), intent(inout) :: this

    if ( allocated(this%temp) ) deallocate( this%temp )
    if ( allocated(this%torr) ) deallocate( this%torr )
    if ( allocated(this%mech) ) deallocate( this%mech )
    
    if ( allocated(this%u_dn) ) deallocate( this%u_dn )
    if ( allocated(this%u_up) ) deallocate( this%u_up )
    if ( allocated(this%u_I2) ) deallocate( this%u_I2 )
    if ( allocated(this%u_C ) ) deallocate( this%u_C  )
    
    if ( allocated(this%t_dn) ) deallocate( this%t_dn )
    if ( allocated(this%t_up) ) deallocate( this%t_up )
    
    if ( allocated(this%v_up) ) deallocate( this%v_up )
    if ( allocated(this%v_dn) ) deallocate( this%v_dn )

  end subroutine deallocate_solution_sub

end module Solution
