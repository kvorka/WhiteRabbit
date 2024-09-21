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
  use Boundaries
  use Matrices
  use Tidal_heating
  implicit none
  
  type, abstract, public :: T_physicalObject
    character(len=6)               :: rheology, scaling
    character(len=5)               :: thermal_bnd, mechanic_bnd, grid_type
    logical                        :: noharm, noobj
    integer                        :: nd, jmax, jms, jms2, jmv, jmt, n_iter, poc
    real(kind=dbl)                 :: t, dt, cf, ab, rd, ru, r_ud, D_ud, gd, gu, Pr, Ra, Ek, St, Cl, Ds, Raf, Ramu, Rad, Rau
    integer,           allocatable :: j_indx(:)
    complex(kind=dbl), allocatable :: flux_up(:), rsph1(:,:), rsph2(:,:), rtorr(:,:), rtemp(:,:), &
                                    & nsph1(:,:), nsph2(:,:), ntorr(:,:), ntemp(:,:)
    
    type(T_radialGrid)   :: rad_grid
    type(T_lateralGrid)  :: lat_grid
    type(T_gravity)      :: gravity
    type(T_Mparams)      :: mparams
    type(T_matrices)     :: mat
    type(T_solution)     :: sol
    type(T_boundaries)   :: bnd
    type(T_tidalHeating) :: tdheat
    
    contains
    
    procedure, pass :: init_objects_sub       => init_objects_sub
    procedure, pass :: deallocate_objects_sub => deallocate_objects_sub
    
    !Time stepping
    procedure, pass :: set_dt_sub
    
    !Material parameters
    procedure, pass :: lambda_r_fn, cp_r_fn, visc_r_fn, alpha_r_fn
    procedure, pass :: lambda_rr_fn, cp_rr_fn, visc_rr_fn, alpha_rr_fn
    
    !Variables :: Thermal solution
    procedure, pass :: htide_r_fn, htide_ir_ijm_sub
    procedure, pass :: htide_rr_fn
    procedure, pass :: temp_r_fn, temp_r_ijm_sub, temp_ir_jm_sub, dT_dr_r_fn, dT_dr_r_ijm_sub, gradT_r_ijml_sub
    procedure, pass :: temp_rr_fn, temp_rr_ijm_sub, temp_irr_jm_sub, dT_dr_rr_fn, dT_dr_rr_ijm_sub, gradT_rr_ijml_sub
    procedure, pass :: q_r_fn, q_r_ijml_sub
    procedure, pass :: qr_r_fn, qr_r_ijm_sub, qr_ir_jm_sub
    procedure, pass :: q_rr_fn, q_rr_ijml_sub, dq_dr_rr_fn, dq_dr_rr_ijml_sub, divq_rr_ijm_sub
    procedure, pass :: qr_rr_fn, qr_rr_ijm_sub, qr_irr_jm_sub
    
    !Variables mechanical part of solution
    procedure, pass :: v_r_fn, v_r_ijml_sub
    procedure, pass :: v_rr_fn, v_rr_ijml_sub
    procedure, pass :: vr_r_fn, vr_rr_fn, vr_r_jm_sub, vr_rr_jm_sub, dv_dr_rr_jml_sub
    
    !Matrices, equations, solvers
    procedure, pass :: init_eq_mech_sub, init_eq_torr_sub, init_eq_temp_sub
    procedure, pass :: mat_temp_fn, mat_mech_fn, mat_torr_fn
    procedure, pass :: prepare_mat_mech_sub, prepare_mat_temp_sub, prepare_mat_torr_sub
    procedure, pass :: solve_temp_sub, solve_torr_sub, solve_mech_sub
    
    !Forces and non-linear terms
    procedure, pass :: volume_heating_fn, tidal_heating_4_sub
    procedure, pass :: global_rotation_sub
    procedure, pass :: coriolis_sub, coriolis_rr_jml_sub
    procedure, pass :: buoy_rr_fn, buoy_rr_jml_sub, er_buoy_rr_jm_sub
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
    
    !Interfaces :: Variables temperature
    module pure complex(kind=dbl) function temp_r_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function temp_r_fn
    
    module pure subroutine temp_r_ijm_sub(this, ir, temp_r_ijm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: temp_r_ijm(:)
    end subroutine temp_r_ijm_sub
    
    module pure subroutine temp_ir_jm_sub(this, ijm, temp_ir_jm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ijm
      complex(kind=dbl),       intent(out) :: temp_ir_jm(:)
    end subroutine temp_ir_jm_sub
    
    module pure complex(kind=dbl) function dT_dr_r_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function dT_dr_r_fn
    
    module pure subroutine dT_dr_r_ijm_sub(this, ir, dT_dr_r)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: dT_dr_r(:)
    end subroutine dT_dr_r_ijm_sub
    
    module pure subroutine gradT_r_ijml_sub(this, ir, gradT, sgn)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir, sgn
      complex(kind=dbl),       intent(out) :: gradT(:)
    end subroutine gradT_r_ijml_sub
    
    module pure complex(kind=dbl) function temp_rr_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function temp_rr_fn
    
    module pure subroutine temp_rr_ijm_sub(this, ir, temp_rr_ijm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: temp_rr_ijm(:)
    end subroutine temp_rr_ijm_sub
    
    module pure subroutine temp_irr_jm_sub(this, ijm, temp_irr_jm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ijm
      complex(kind=dbl),       intent(out) :: temp_irr_jm(:)
    end subroutine temp_irr_jm_sub
    
    module pure complex(kind=dbl) function dT_dr_rr_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function dT_dr_rr_fn
    
    module pure subroutine dT_dr_rr_ijm_sub(this, ir, T, dT)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: T(:), dT(:)
    end subroutine dT_dr_rr_ijm_sub
    
    module pure subroutine gradT_rr_ijml_sub(this, ir, T, gradT, sgn)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir, sgn
      complex(kind=dbl),       intent(out) :: T(:), gradT(:)
    end subroutine gradT_rr_ijml_sub
    
    !! Interfaces :: variables heat flux
    module pure complex(kind=dbl) function q_r_fn(this, ir, il, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, il, ijm
    end function q_r_fn
    
    module pure subroutine q_r_ijml_sub(this, ir, q_r_ijml)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: q_r_ijml(:)
    end subroutine q_r_ijml_sub
    
    module pure complex(kind=dbl) function qr_r_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function qr_r_fn
    
    module pure subroutine qr_r_ijm_sub(this, ir, qr_r_ijm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: qr_r_ijm(*)
    end subroutine qr_r_ijm_sub
    
    module pure subroutine qr_ir_jm_sub(this, ijm, qr_ir_jm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ijm
      complex(kind=dbl),       intent(out) :: qr_ir_jm(*)
    end subroutine qr_ir_jm_sub
    
    module pure complex(kind=dbl) function q_rr_fn(this, ir, il, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, il, ijm
    end function q_rr_fn
    
    module pure subroutine q_rr_ijml_sub(this, ir, q_rr_ijml)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: q_rr_ijml(:)
    end subroutine q_rr_ijml_sub
    
    module pure complex(kind=dbl) function dq_dr_rr_fn(this, ir, il, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, il, ijm
    end function dq_dr_rr_fn
    
    module pure subroutine dq_dr_rr_ijml_sub(this, ir, dq_dr_rr_ijml)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: dq_dr_rr_ijml(:)
    end subroutine dq_dr_rr_ijml_sub
    
    module pure subroutine divq_rr_ijm_sub(this, ir, divq_rr_ijm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: divq_rr_ijm(:)
    end subroutine divq_rr_ijm_sub
    
    module pure complex(kind=dbl) function qr_rr_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function qr_rr_fn
    
    module pure subroutine qr_rr_ijm_sub(this, ir, qr_rr_ijm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: qr_rr_ijm(*)
    end subroutine qr_rr_ijm_sub
    
    module pure subroutine qr_irr_jm_sub(this, ijm, qr_irr_jm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ijm
      complex(kind=dbl),       intent(out) :: qr_irr_jm(2:this%nd)
    end subroutine qr_irr_jm_sub
    
    !Interfaces :: Variables velocity
    module pure complex(kind=dbl) function v_r_fn(this, ir, il, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, il, ijm
    end function v_r_fn
    
    module pure subroutine v_r_ijml_sub(this, ir, v_r_ijml)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: v_r_ijml(:)
    end subroutine v_r_ijml_sub
    
    module pure complex(kind=dbl) function v_rr_fn(this, ir, il, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, il, ijm
    end function v_rr_fn
    
    module pure subroutine v_rr_ijml_sub(this, ir, v_rr_ijml)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: v_rr_ijml(:)
    end subroutine v_rr_ijml_sub
    
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
      complex(kind=dbl),       intent(out) :: vr_jm(:)
    end subroutine vr_r_jm_sub
    
    module pure subroutine vr_rr_jm_sub(this, ir, vr_jm)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: vr_jm(:)
    end subroutine vr_rr_jm_sub
    
    module pure subroutine dv_dr_rr_jml_sub(this, ir, v, dv)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: dv(:), v(:)
    end subroutine dv_dr_rr_jml_sub
    
    !Interfaces :: material parameters
    module pure real(kind=dbl) function lambda_r_fn(this, ir)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir
    end function lambda_r_fn
    
    module pure real(kind=dbl) function cp_r_fn(this, ir)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir
    end function cp_r_fn
    
    module pure real(kind=dbl) function visc_r_fn(this, ir)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir
    end function visc_r_fn
    
    module pure real(kind=dbl) function alpha_r_fn(this, ir)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir
    end function alpha_r_fn
    
    module pure real(kind=dbl) function lambda_rr_fn(this, ir)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir
    end function lambda_rr_fn
    
    module pure real(kind=dbl) function cp_rr_fn(this, ir)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir
    end function cp_rr_fn
    
    module pure real(kind=dbl) function visc_rr_fn(this, ir)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir
    end function visc_rr_fn
    
    module pure real(kind=dbl) function alpha_rr_fn(this, ir)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir
    end function alpha_rr_fn
    
    !Interfaces :: tidal heating
    module pure complex(kind=dbl) function htide_r_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function htide_r_fn
    
    module pure subroutine htide_ir_ijm_sub(this, htide)
      class(T_physicalObject), intent(in)    :: this
      complex(kind=dbl),       intent(inout) :: htide(:,:)
    end subroutine htide_ir_ijm_sub
    
    module pure complex(kind=dbl) function htide_rr_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function htide_rr_fn
    
    !Interfaces :: output
    module subroutine vypis_sub(this, filenum, path, quantity)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: filenum
      character(len=*),        intent(in) :: path, quantity
    end subroutine vypis_sub
    
    !Interfaces :: to be continued
    module pure subroutine coriolis_rr_jml_sub(this, v, coriolis)
      class(T_physicalObject), intent(in)    :: this
      complex(kind=dbl),       intent(in)    :: v(:)
      complex(kind=dbl),       intent(inout) :: coriolis(:,:)
    end subroutine coriolis_rr_jml_sub
    
    module pure complex(kind=dbl) function buoy_rr_fn(this, ir, il, ijm, sgn)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, il, ijm, sgn
    end function buoy_rr_fn
    
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
    
    module subroutine tidal_heating_4_sub(this)
      class(T_physicalObject), intent(inout) :: this
    end subroutine tidal_heating_4_sub
    
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
    
    module real(kind=dbl) function reynolds_fn(this, choice)
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