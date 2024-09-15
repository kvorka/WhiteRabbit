module PhysicalObject
  use Math
  use sph_indexing
  use sph_unitvec_op
  use sph_norms
  use SphericalHarmonics
  use RadialGrid
  use Gravity
  use Mparams
  use Solution
  use Matrices
  implicit none
  
  type, abstract, public :: T_physicalObject
    character(len=6)               :: rheology, scaling
    character(len=5)               :: thermal_bnd, mechanic_bnd, grid_type
    logical                        :: noharm, noobj
    integer                        :: nd, jmax, jms, jms2, jmv, jmt, n_iter, poc
    real(kind=dbl)                 :: t, dt, cf, ab, rd, ru, r_ud, D_ud, gd, gu, Pr, Ra, Ek, St, Cl, Ds, Raf, Ramu, Rad, Rau
    integer,           allocatable :: j_indx(:)
    complex(kind=dbl), allocatable :: flux_up(:), htide(:,:), rsph1(:,:), rsph2(:,:), rtorr(:,:), rtemp(:,:), &
                                    &                         nsph1(:,:), nsph2(:,:), ntorr(:,:), ntemp(:,:)
    
    type(T_radialGrid)  :: rad_grid
    type(T_lateralGrid) :: lat_grid
    type(T_gravity)     :: gravity
    type(T_Mparams)     :: mparams
    type(T_matrices)    :: mat
    type(T_solution)    :: sol
    
    contains
    
    procedure, pass :: init_objects_sub       => init_objects_sub
    procedure, pass :: deallocate_objects_sub => deallocate_objects_sub
    
    !Time stepping
    procedure, pass :: set_dt_sub
    
    !Material parameters
    procedure, pass :: lambda_fn, cp_fn, visc_fn, alpha_fn
    
    !Variables
    procedure, pass :: htide_fn
    procedure, pass :: temp_r_fn, temp_rr_fn, qr_r_fn, mgradT_rr_jml_sub
    procedure, pass :: vr_r_fn, vr_rr_fn, vr_r_jm_sub, vr_rr_jm_sub, dv_dr_rr_jml_sub
    
    !Matrices, equations, solvers
    procedure, pass :: init_eq_mech_sub, init_eq_torr_sub, init_eq_temp_sub
    procedure, pass :: mat_temp_fn, mat_mech_fn, mat_torr_fn
    procedure, pass :: prepare_mat_mech_sub, prepare_mat_temp_sub, prepare_mat_torr_sub
    procedure, pass :: solve_temp_sub, solve_torr_sub, solve_mech_sub
    
    !Forces and non-linear terms
    procedure, pass :: volume_heating_fn
    procedure, pass :: global_rotation_sub
    procedure, pass :: coriolis_sub, coriolis_rr_jml_sub
    procedure, pass :: buoy_rr_jml_sub, er_buoy_rr_jm_sub
    procedure, pass :: viscdissip_power_fn, buoyancy_power_fn, bottombnd_power_fn, upperbnd_power_fn
    procedure, pass :: coriolis_vgradv_sub, mvgradT_sub, fullnl_sub
    
    !Output, control measures
    procedure, pass :: vypis_sub
    procedure, pass :: reynolds_fn, nuss_fn
    procedure, pass :: laws_temp_fn, laws_mech_fn, laws_force_fn
    
  end type T_physicalObject
  
  interface
    module subroutine init_objects_sub(this, nd, jmax, r_ud, rgrid, gmod, g, noobj, noharm)
      class(T_physicalObject),    intent(inout) :: this
      integer,                    intent(in)    :: nd, jmax
      real(kind=dbl),             intent(in)    :: r_ud, g
      character(len=*),           intent(in)    :: rgrid, gmod
      logical,          optional, intent(in)    :: noobj, noharm
    end subroutine init_objects_sub
    
    module subroutine set_dt_sub(this)
      class(T_physicalObject), intent(inout) :: this
    end subroutine set_dt_sub
    
    module subroutine deallocate_objects_sub(this)
      class(T_physicalObject), intent(inout) :: this
    end subroutine deallocate_objects_sub
    
    module subroutine vypis_sub(this, filenum, path, quantity)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: filenum
      character(len=*),        intent(in) :: path, quantity
    end subroutine vypis_sub
    
    module pure real(kind=dbl) function lambda_fn(this, i)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: i
    end function lambda_fn
    
    module pure real(kind=dbl) function cp_fn(this, i)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: i
    end function cp_fn
    
    module pure real(kind=dbl) function visc_fn(this, i)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: i
    end function visc_fn
    
    module pure real(kind=dbl) function alpha_fn(this, i)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: i
    end function alpha_fn
    
    module pure complex(kind=dbl) function htide_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function htide_fn
    
    module pure complex(kind=dbl) function temp_r_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function temp_r_fn
    
    module pure complex(kind=dbl) function temp_rr_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function temp_rr_fn
    
    module pure complex(kind=dbl) function vr_r_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function vr_r_fn
    
    module pure complex(kind=dbl) function vr_rr_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function vr_rr_fn
    
    module pure subroutine vr_r_jm_sub(this, ir, vr_jm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: vr_jm(*)
    end subroutine vr_r_jm_sub
    
    module pure subroutine vr_rr_jm_sub(this, ir, vr_jm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: vr_jm(*)
    end subroutine vr_rr_jm_sub
    
    module pure complex(kind=dbl) function qr_r_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function qr_r_fn
    
    module pure subroutine dv_dr_rr_jml_sub(this, ir, v, dv)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: dv(:), v(:)
    end subroutine dv_dr_rr_jml_sub
    
    module pure subroutine mgradT_rr_jml_sub(this, ir, T, gradT)
      class(T_physicalObject),     intent(in)  :: this
      integer,                     intent(in)  :: ir
      complex(kind=dbl), optional, intent(out) :: T(:)
      complex(kind=dbl),           intent(out) :: gradT(:)
    end subroutine mgradT_rr_jml_sub
    
    module pure subroutine coriolis_rr_jml_sub(this, v, coriolis)
      class(T_physicalObject), intent(in)    :: this
      complex(kind=dbl),       intent(in)    :: v(:)
      complex(kind=dbl),       intent(inout) :: coriolis(:,:)
    end subroutine coriolis_rr_jml_sub
    
    module pure subroutine er_buoy_rr_jm_sub(this, ir, force)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: force(*)
    end subroutine er_buoy_rr_jm_sub
    
    module pure subroutine buoy_rr_jml_sub(this, ir, T, force)
      class(T_physicalObject), intent(in)    :: this
      integer,                 intent(in)    :: ir
      complex(kind=dbl),       intent(in)    :: T(:)
      complex(kind=dbl),       intent(inout) :: force(:,:)
    end subroutine buoy_rr_jml_sub
    
    module pure subroutine global_rotation_sub(this)
      class(T_physicalObject), intent(inout) :: this
    end subroutine global_rotation_sub
    
    module subroutine mvgradT_sub(this, i, mvgradT)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: i
      complex(kind=dbl),       intent(out) :: mvgradT(:)
    end subroutine mvgradT_sub
    
    module subroutine fullnl_sub(this, i)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: i
    end subroutine fullnl_sub

    module subroutine coriolis_vgradv_sub(this, i)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: i
    end subroutine coriolis_vgradv_sub

    module pure subroutine coriolis_sub(this, i)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: i
    end subroutine coriolis_sub
    
    module pure function viscdissip_power_fn(this) result(power)
      class(T_physicalObject), intent(in)  :: this
      real(kind=dbl)                       :: power
    end function viscdissip_power_fn
    
    module pure function buoyancy_power_fn(this) result(power)
      class(T_physicalObject), intent(in)  :: this
      real(kind=dbl)                       :: power
    end function buoyancy_power_fn
    
    module pure function bottombnd_power_fn(this) result(power)
      class(T_physicalObject), intent(in) :: this
      real(kind=dbl)                      :: power
    end function bottombnd_power_fn
    
    module pure function upperbnd_power_fn(this) result(power)
      class(T_physicalObject), intent(in) :: this
      real(kind=dbl)                      :: power
    end function upperbnd_power_fn
    
    module pure function mat_temp_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),        allocatable  :: matica(:,:)
    end function mat_temp_fn
    
    module pure function mat_torr_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),        allocatable  :: matica(:,:)
    end function mat_torr_fn
    
    module pure function mat_mech_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),        allocatable  :: matica(:,:)
    end function mat_mech_fn
    
    module pure function matica_mech_chb_viscos_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),         allocatable :: matica(:,:)
    end function matica_mech_chb_viscos_fn
    
    module pure function matica_mech_chb_christ_viscos_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),         allocatable :: matica(:,:)
    end function matica_mech_chb_christ_viscos_fn
    
    module pure function matica_mech_hom_viscos_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),         allocatable :: matica(:,:)
    end function matica_mech_hom_viscos_fn
    
    module pure function matica_mech_chb_viscel_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),         allocatable :: matica(:,:)
    end function matica_mech_chb_viscel_fn
    
    module pure function matica_mech_hom_viscel_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),         allocatable :: matica(:,:)
    end function matica_mech_hom_viscel_fn
    
    module pure function matica_temp_hom_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),         allocatable :: matica(:,:)
    end function matica_temp_hom_fn
    
    module pure function matica_temp_chb_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),         allocatable :: matica(:,:)
    end function matica_temp_chb_fn
    
    module pure function matica_torr_chb_viscos_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),         allocatable :: matica(:,:)
    end function matica_torr_chb_viscos_fn
    
    module pure function matica_torr_chb_christ_viscos_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),         allocatable :: matica(:,:)
    end function matica_torr_chb_christ_viscos_fn
    
    module subroutine init_eq_temp_sub(this,rhs,nl)
      class(T_physicalObject), intent(inout) :: this
      logical,                 intent(in)    :: rhs, nl
    end subroutine init_eq_temp_sub
    
    module subroutine init_eq_torr_sub(this,rhs,nl)
      class(T_physicalObject), intent(inout) :: this
      logical,                 intent(in)    :: rhs, nl
    end subroutine init_eq_torr_sub
    
    module subroutine init_eq_mech_sub(this,rhs,nl)
      class(T_physicalObject), intent(inout) :: this
      logical,                 intent(in)    :: rhs, nl
    end subroutine init_eq_mech_sub
    
    module subroutine prepare_mat_temp_sub(this, ijstart, ijend)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: ijstart, ijend
    end subroutine prepare_mat_temp_sub
    
    module subroutine prepare_mat_torr_sub(this, ijstart, ijend)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: ijstart, ijend
    end subroutine prepare_mat_torr_sub
    
    module subroutine prepare_mat_mech_sub(this, ijstart, ijend)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: ijstart, ijend
    end subroutine prepare_mat_mech_sub
    
    module subroutine solve_temp_sub(this, ijmstart, ijmend, ijmstep, rematrix, matxsol)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: ijmstart, ijmend, ijmstep
      logical,                 intent(in)    :: rematrix, matxsol
    end subroutine solve_temp_sub
    
    module subroutine solve_torr_sub(this, ijmstart, ijmend, ijmstep, rematrix, matxsol)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: ijmstart, ijmend, ijmstep
      logical,                 intent(in)    :: rematrix, matxsol
    end subroutine solve_torr_sub
    
    module subroutine solve_mech_sub(this, ijmstart, ijmend, ijmstep, rematrix, matxsol)
      class(T_physicalObject), intent(inout)        :: this
      integer,                 intent(in)           :: ijmstart, ijmend, ijmstep
      logical,                 intent(in)           :: rematrix, matxsol
    end subroutine solve_mech_sub
    
    module pure real(kind=dbl) function nuss_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function nuss_fn
    
    module pure real(kind=dbl) function reynolds_fn(this, choice)
      class(T_physicalObject), intent(in)           :: this
      character(len=*),        intent(in), optional :: choice
    end function reynolds_fn
    
    module pure real(kind=dbl) function volume_heating_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function volume_heating_fn
    
    module real(kind=dbl) function laws_mech_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function laws_mech_fn
    
    module real(kind=dbl) function laws_temp_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function laws_temp_fn
    
    module real(kind=dbl) function laws_force_fn(this, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ijm
    end function laws_force_fn
  end interface

end module PhysicalObject