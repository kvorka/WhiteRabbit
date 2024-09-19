module Solution
  use Math
  use Nulify
  implicit none
  
  type, public :: T_solution
    logical                        :: inittemp, initsfer, inittorr
    integer                        :: nd, jmax, jms, jmv, jmt
    complex(kind=dbl), allocatable :: u_dn(:), u_up(:), u_I2(:), u_C(:), v_dn(:), v_up(:), t_dn(:), t_up(:), &
                                    & mech(:,:), temp(:,:), torr(:,:)
    
    contains
    
    procedure :: init_sub       => init_solution_sub
    procedure :: nulify_sub     => nulify_solution_sub
    procedure :: deallocate_sub => deallocate_solution_sub
    
    procedure :: init_stemp_sub, temp_fn, flux_fn
    procedure :: temp_rr_many1_sub
    procedure :: temp_jm_many1_sub, temp_jm_many2_sub, temp_jm_many3_sub, temp_jm_many4_sub
    procedure :: flux_r_many1_sub
    procedure :: flux_jml_many1_sub, flux_jml_many2_sub, flux_jml_many3_sub, flux_jml_many4_sub
    
    procedure :: init_storr_sub, init_smech_sub, velocity_fn, deviatoric_stress_fn
    procedure :: velocity_rr_many1_sub
    procedure :: conv_velocity_jml_sub, velocity_jml_many1_sub, velocity_jml_many2_sub, velocity_jml_many3_sub
    procedure :: deviatoric_stress_jml2_fn
    
    procedure :: init_layers_sub, init_layer_u_sub
    
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
    
    !! Interfaces :: temperature
    module pure subroutine init_stemp_sub(this)
      class(T_solution), intent(inout) :: this
    end subroutine init_stemp_sub
    
    module pure complex(kind=dbl) function temp_fn(this, ir, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, ijm
    end function temp_fn
    
    module pure subroutine temp_rr_many1_sub(this, ijm, temp1)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ijm
      complex(kind=dbl), intent(out) :: temp1(:)
    end subroutine temp_rr_many1_sub
    
    module pure subroutine temp_jm_many1_sub(this, ir, temp1)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: temp1(:)
    end subroutine temp_jm_many1_sub
    
    module pure subroutine temp_jm_many2_sub(this, ir, temp1, temp2)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: temp1(:), temp2(:)
    end subroutine temp_jm_many2_sub
    
    module pure subroutine temp_jm_many3_sub(this, ir, temp1, temp2, temp3)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: temp1(:), temp2(:), temp3(:)
    end subroutine temp_jm_many3_sub
    
    module pure subroutine temp_jm_many4_sub(this, ir, temp1, temp2, temp3, temp4)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: temp1(:), temp2(:), temp3(:), temp4(:)
    end subroutine temp_jm_many4_sub
    
    module pure complex(kind=dbl) function flux_fn(this, ir, il, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, il, ijm
    end function flux_fn
    
    module pure subroutine flux_r_many1_sub(this, ijm, flux1)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ijm
      complex(kind=dbl), intent(out) :: flux1(:,:)
    end subroutine flux_r_many1_sub
    
    module pure subroutine flux_jml_many1_sub(this, ir, flux1)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: flux1(:)
    end subroutine flux_jml_many1_sub
    
    module pure subroutine flux_jml_many2_sub(this, ir, flux1, flux2)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: flux1(:), flux2(:)
    end subroutine flux_jml_many2_sub
    
    module pure subroutine flux_jml_many3_sub(this, ir, flux1, flux2, flux3)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: flux1(:), flux2(:), flux3(:)
    end subroutine flux_jml_many3_sub
    
    module pure subroutine flux_jml_many4_sub(this, ir, flux1, flux2, flux3, flux4)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: flux1(:), flux2(:), flux3(:), flux4(:)
    end subroutine flux_jml_many4_sub
    
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
    
    module pure complex(kind=dbl) function velocity_fn(this, ir, il, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, il, ijm
    end function velocity_fn
    
    module pure subroutine velocity_rr_many1_sub(this, ijm, velocity1)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ijm
      complex(kind=dbl), intent(out) :: velocity1(:,:)
    end subroutine velocity_rr_many1_sub
    
    module pure subroutine conv_velocity_jml_sub(this, ir, velocity)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: velocity(:)
    end subroutine conv_velocity_jml_sub
    
    module pure subroutine velocity_jml_many1_sub(this, ir, velocity1)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: velocity1(:)
    end subroutine velocity_jml_many1_sub
    
    module pure subroutine velocity_jml_many2_sub(this, ir, velocity1, velocity2)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: velocity1(:), velocity2(:)
    end subroutine velocity_jml_many2_sub
    
    module pure subroutine velocity_jml_many3_sub(this, ir, velocity1, velocity2, velocity3)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), intent(out) :: velocity1(:), velocity2(:), velocity3(:)
    end subroutine velocity_jml_many3_sub
    
    module pure complex(kind=dbl) function deviatoric_stress_fn(this, ir, il, ijm)
      class(T_solution), intent(in) :: this
      integer,           intent(in) :: ir, il, ijm
    end function deviatoric_stress_fn
    
    module pure function deviatoric_stress_jml2_fn(this, ir) result(dstress)
      class(T_solution), intent(in)  :: this
      integer,           intent(in)  :: ir
      complex(kind=dbl), allocatable :: dstress(:)
    end function deviatoric_stress_jml2_fn
  end interface
  
end module Solution