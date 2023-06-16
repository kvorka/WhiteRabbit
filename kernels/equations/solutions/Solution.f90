module Solution
  use Math
  implicit none

  type, public :: T_solution 
    integer                        :: nd, jmax, jms, jmv, jmt
    complex(kind=dbl), allocatable :: u_dn(:), u_up(:), u_I2(:), u_C(:), v_dn(:), v_up(:), t_dn(:), t_up(:), &
                                    & mech(:,:), temp(:,:), torr(:,:)
    
    contains
    
    procedure, pass :: init_sub       => init_solution_sub
    procedure, pass :: deallocate_sub => deallocate_solution_sub
    
    procedure, pass :: init_stemp_sub
    procedure, pass :: init_storr_sub
    procedure, pass :: init_smech_sub
    procedure, pass :: init_layers_sub
    procedure, pass :: init_layer_u_sub
    procedure, pass :: temp_fn
    procedure, pass :: temp2_fn
    procedure, pass :: flux_fn
    procedure, pass :: velocity_fn
    procedure, pass :: deviatoric_stress_fn
    procedure, pass :: temp_i_fn
    procedure, pass :: flux_i_fn
    procedure, pass :: velocity_i_fn
    procedure, pass :: deviatoric_stress_i_fn
    procedure, pass :: temp_jm_fn
    procedure, pass :: temp_jm_sub
    procedure, pass :: flux_jml_fn
    procedure, pass :: flux_jml_sub
    procedure, pass :: flux_jml_many_sub
    procedure, pass :: velocity_jml_fn
    procedure, pass :: velocity_jml_sub
    procedure, pass :: velocity_jml_many_sub
    procedure, pass :: conv_velocity_jml_fn
    procedure, pass :: deviatoric_stress_jml2_fn
    
  end type T_solution

  interface
    module subroutine init_solution_sub(this, nd, jmax)
      class(T_solution), intent(inout) :: this
      integer,           intent(in)    :: nd, jmax
    end subroutine init_solution_sub
    
    module subroutine deallocate_solution_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine deallocate_solution_sub
    
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
    
    module pure complex(kind=dbl) function temp_fn(this, ir, ij, im)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, ij, im
    end function temp_fn
    
    module pure complex(kind=dbl) function temp2_fn(this, ir, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, ijm
    end function temp2_fn
    
    module pure complex(kind=dbl) function flux_fn(this, ir, ij, im, il)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, ij, im, il
    end function flux_fn

    module pure complex(kind=dbl) function flux2_fn(this, ir, il, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, ijm, il
    end function flux2_fn
    
    module pure complex(kind=dbl) function velocity_fn(this, ir, ij, im, il)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, ij, im, il
    end function velocity_fn

    module pure complex(kind=dbl) function velocity2_fn(this, ir, il, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, ijm, il
    end function velocity2_fn
    
    module pure complex(kind=dbl) function deviatoric_stress_fn(this, ir, ij, im, il)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, ij, im, il
    end function deviatoric_stress_fn

    module pure complex(kind=dbl) function deviatoric_stress2_fn(this, ir, il, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, ijm, il
    end function deviatoric_stress2_fn
    
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

    module pure subroutine flux_jml_many_sub(this, i, temp2, flux1, flux2)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: i
      complex(kind=dbl), intent(out) :: temp2(:), flux1(:), flux2(:)
    end subroutine flux_jml_many_sub
    
    module pure function velocity_jml_fn(this, i) result(velocity)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: i
      complex(kind=dbl)             :: velocity(this%jmv)
    end function velocity_jml_fn

    module pure subroutine velocity_jml_many_sub(this, i, velocity1, velocity2, velocity3)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: i
      complex(kind=dbl), intent(out) :: velocity1(:), velocity2(:), velocity3(:)
    end subroutine velocity_jml_many_sub
    
    module pure subroutine temp_jm_sub(this, i, temp)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: i
      complex(kind=dbl), intent(out) :: temp(:)
    end subroutine temp_jm_sub
    
    module pure subroutine flux_jml_sub(this, i, flux)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: i
      complex(kind=dbl), intent(out) :: flux(:)
    end subroutine flux_jml_sub
    
    module pure subroutine velocity_jml_sub(this, i, velocity)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: i
      complex(kind=dbl), intent(out) :: velocity(:)
    end subroutine velocity_jml_sub
    
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
  
end module Solution
