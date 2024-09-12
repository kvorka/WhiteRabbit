module Solution
  use Math
  implicit none
  
  type, public :: T_solution
    logical                        :: inittemp, initsfer, inittorr
    integer                        :: nd, jmax, jms, jmv, jmt
    complex(kind=dbl), allocatable :: u_dn(:), u_up(:), u_I2(:), u_C(:), v_dn(:), v_up(:), t_dn(:), t_up(:), &
                                    & mech(:,:), temp(:,:), torr(:,:), visc(:,:)
    
    contains
    
    procedure :: init_sub       => init_solution_sub
    procedure :: nulify_sub     => nulify_solution_sub
    procedure :: deallocate_sub => deallocate_solution_sub
    procedure :: init_stemp_sub, init_storr_sub, init_smech_sub, init_layers_sub, init_layer_u_sub, temp_fn, flux_fn,        &
               & velocity_fn, deviatoric_stress_fn, temp_i_fn, flux_i_fn, velocity_i_fn, deviatoric_stress_i_fn, temp_jm_fn, &
               & temp_jm_sub, flux_jml_fn, flux_jml_sub, flux_jml_many_sub, velocity_jml_fn, velocity_jml_sub,               &
               & velocity_jml_many_sub, conv_velocity_jml_fn, deviatoric_stress_jml2_fn, init_visc_sub
    
  end type T_solution
  
  interface
    module pure subroutine init_solution_sub(this, nd, jmax)
      class(T_solution), intent(inout) :: this
      integer,           intent(in)    :: nd, jmax
    end subroutine init_solution_sub
    
    module pure subroutine nulify_solution_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine nulify_solution_sub
    
    module pure subroutine deallocate_solution_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine deallocate_solution_sub
    
    module pure subroutine init_stemp_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_stemp_sub
    
    module pure subroutine init_storr_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_storr_sub
    
    module pure subroutine init_smech_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_smech_sub
    
    module pure subroutine init_layers_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_layers_sub
    
    module pure subroutine init_layer_u_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_layer_u_sub
    
    module pure subroutine init_visc_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_visc_sub
    
    module pure complex(kind=dbl) function temp_fn(this, ir, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, ijm
    end function temp_fn
    
    module pure complex(kind=dbl) function flux_fn(this, ir, il, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, il, ijm
    end function flux_fn
    
    module pure complex(kind=dbl) function velocity_fn(this, ir, il, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, il, ijm
    end function velocity_fn
    
    module pure complex(kind=dbl) function deviatoric_stress_fn(this, ir, il, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, il, ijm
    end function deviatoric_stress_fn
    
    module pure function temp_i_fn(this, ijm) result(temp)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ijm
      complex(kind=dbl), allocatable :: temp(:)
    end function temp_i_fn
    
    module pure function flux_i_fn(this, il, ijm) result(flux)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: il, ijm
      complex(kind=dbl), allocatable :: flux(:)
    end function flux_i_fn
    
    module pure function velocity_i_fn(this, il, ijm) result(velocity)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: il, ijm
      complex(kind=dbl), allocatable :: velocity(:)
    end function velocity_i_fn
    
    module pure function deviatoric_stress_i_fn(this, il, ijm) result(dstress)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: il, ijm
      complex(kind=dbl), allocatable :: dstress(:)
    end function deviatoric_stress_i_fn
    
    module pure function temp_jm_fn(this, ir) result(temp)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), allocatable :: temp(:)
    end function temp_jm_fn
    
    module pure function flux_jml_fn(this, ir) result(flux)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), allocatable :: flux(:)
    end function flux_jml_fn
    
    module pure subroutine flux_jml_many_sub(this, ir, temp2, flux1, flux2)
      class(T_solution),           intent(in)  :: this
      integer,                     intent(in)  :: ir
      complex(kind=dbl), optional, intent(out) :: temp2(:)
      complex(kind=dbl),           intent(out) :: flux1(:), flux2(:)
    end subroutine flux_jml_many_sub
    
    module pure function velocity_jml_fn(this, ir) result(velocity)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), allocatable :: velocity(:)
    end function velocity_jml_fn
    
    module pure subroutine velocity_jml_many_sub(this, ir, velocity1, velocity2, velocity3)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: velocity1(:), velocity2(:), velocity3(:)
    end subroutine velocity_jml_many_sub
    
    module pure subroutine temp_jm_sub(this, ir, temp)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: temp(:)
    end subroutine temp_jm_sub
    
    module pure subroutine flux_jml_sub(this, ir, flux)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: flux(:)
    end subroutine flux_jml_sub
    
    module pure subroutine velocity_jml_sub(this, ir, velocity)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: velocity(:)
    end subroutine velocity_jml_sub
    
    module pure function conv_velocity_jml_fn(this, ir) result(velocity)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), allocatable :: velocity(:)
    end function conv_velocity_jml_fn
    
    module pure function deviatoric_stress_jml2_fn(this, ir) result(dstress)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), allocatable :: dstress(:)
    end function deviatoric_stress_jml2_fn
  end interface
  
end module Solution