module PhysicalObject
  use Math
  use Spherical_func
  use SphericalHarmonics
  use RadialGrid
  use Gravity
  use Solution
  use Matrices
  implicit none
  
  type, abstract, public :: T_physicalObject
    character(len=6)               :: rheology
    character(len=5)               :: thermal_bnd
    character(len=5)               :: mechanic_bnd
    character(len=5)               :: grid_type
    character(len=6)               :: scaling
    logical                        :: noharm, noobj
    integer                        :: nd, jmax, jms, jmv
    integer                        :: n_iter, poc
    real(kind=dbl)                 :: t, dt
    real(kind=dbl)                 :: rd, ru, r_ud, D_ud
    real(kind=dbl)                 :: Pr, Ra, Ek, Ds, Raf, Ramu, Rad, Rau
    integer,           allocatable :: j_indx(:)
    complex(kind=dbl), allocatable :: flux_up(:), htide(:,:)
    complex(kind=dbl), allocatable :: rsph1(:,:), rsph2(:,:), rtorr(:,:), rtemp(:,:) 
    complex(kind=dbl), allocatable :: nsph1(:,:), nsph2(:,:), ntorr(:,:), ntemp(:,:)
    
    type(T_radialGrid)  :: rad_grid
    type(T_lateralGrid) :: lat_grid
    type(T_gravity)     :: gravity
    type(T_matrices)    :: mat
    type(T_solution)    :: sol
    
    contains
    
    procedure, pass :: init_objects_sub       => init_objects_sub
    procedure, pass :: vypis_sub
    procedure, pass :: deallocate_objects_sub => deallocate_objects_sub

    procedure, pass :: lambda_fn
    procedure, pass :: cp_fn
    procedure, pass :: visc_fn
    procedure, pass :: alpha_fn

    procedure, pass :: htide_fn
    procedure, pass :: qr_jm_fn
    procedure, pass :: vr_fn
    procedure, pass :: vr_jm_fn
    procedure, pass :: dv_dr_rrjml_fn
    procedure, pass :: mgradT_rrjml_fn
    procedure, pass :: dv_dr_rr_jml_sub
    procedure, pass :: mgradT_rr_jml_sub

    procedure, pass :: buoy_rr_fn
    procedure, pass :: buoy_rr_jml_fn
    procedure, pass :: coriolis_rr_jml_fn
    procedure, pass :: coriolis_rr_jml_sub
    procedure, pass :: buoy_rr_jml_sub
    procedure, pass :: global_rotation_sub

    procedure, pass :: vgradT_fn
    procedure, pass :: vgradv_fn
    procedure, pass :: fullnl_sub

    procedure, pass :: matica_temp_fn
    procedure, pass :: matica_mech_fn
    procedure, pass :: matica_torr_fn

    procedure, pass :: init_eq_temp_sub
    procedure, pass :: init_eq_mech_sub
    procedure, pass :: init_eq_torr_sub

    procedure, pass :: solve_temp_sub
    procedure, pass :: solve_torr_sub
    procedure, pass :: solve_mech_sub

    procedure, pass :: nuss_fn
    procedure, pass :: reynolds_fn
    procedure, pass :: nonzon_reynolds_fn
    procedure, pass :: volume_heating_fn

    procedure, pass :: laws_mech_fn
    procedure, pass :: laws_temp_fn
    procedure, pass :: laws_force_fn
    
  end type T_physicalObject

  interface
    module subroutine init_objects_sub( this, nd, jmax, r_ud, rgrid, noobj, noharm )
      class(T_physicalObject),    intent(inout) :: this
      integer,                    intent(in)    :: nd, jmax
      real(kind=dbl),             intent(in)    :: r_ud
      character(len=*),           intent(in)    :: rgrid
      logical,          optional, intent(in)    :: noobj, noharm
    end subroutine init_objects_sub
    
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

    module pure complex(kind=dbl) function vr_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function vr_fn

    module pure function qr_jm_fn(this, ir) result(qr)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       allocatable :: qr(:)
    end function qr_jm_fn

    module pure function vr_jm_fn(this, ir) result(vr)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       allocatable :: vr(:)
    end function vr_jm_fn

    module pure function dv_dr_rrjml_fn(this, ir, v) result(dv)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(in)  :: v(:)
      complex(kind=dbl),       allocatable :: dv(:)
    end function dv_dr_rrjml_fn

    module pure function mgradT_rrjml_fn(this, ir) result(gradT)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       allocatable :: gradT(:)
    end function mgradT_rrjml_fn

    module subroutine dv_dr_rr_jml_sub(this, ir, v, dv)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: dv(:), v(:)
    end subroutine dv_dr_rr_jml_sub
    
    module subroutine mgradT_rr_jml_sub(this, ir, T, gradT)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       intent(out) :: T(:), gradT(:)
    end subroutine mgradT_rr_jml_sub
    
    module pure function buoy_rr_fn(this, ir, ijm) result(buoy)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
      complex(kind=dbl)                   :: buoy
    end function buoy_rr_fn
    
    module pure function buoy_rr_jml_fn(this, ir) result(gdrho)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: ir
      complex(kind=dbl),       allocatable :: gdrho(:)
    end function buoy_rr_jml_fn

    module pure function coriolis_rr_jml_fn(this, ir, v) result(coriolis)
      class(T_physicalObject),     intent(in)  :: this
      integer,           optional, intent(in)  :: ir
      complex(kind=dbl), optional, intent(in)  :: v(:)
      complex(kind=dbl),           allocatable :: coriolis(:)
    end function coriolis_rr_jml_fn

    module subroutine coriolis_rr_jml_sub(this, v, coriolis)
      class(T_physicalObject), intent(in)    :: this
      complex(kind=dbl),       intent(in)    :: v(:)
      complex(kind=dbl),       intent(inout) :: coriolis(:,:)
    end subroutine coriolis_rr_jml_sub

    module subroutine buoy_rr_jml_sub(this, ir, T, force)
      class(T_physicalObject), intent(in)    :: this
      integer,                 intent(in)    :: ir
      complex(kind=dbl),       intent(in)    :: T(:)
      complex(kind=dbl),       intent(inout) :: force(:,:)
    end subroutine buoy_rr_jml_sub

    module subroutine global_rotation_sub(this)
      class(T_physicalObject), intent(inout) :: this
    end subroutine global_rotation_sub
    
    module function vgradT_fn(this, i) result(vgradT)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: i
      complex(kind=dbl)                   :: vgradT(this%jms)
    end function vgradT_fn

    module function vgradv_fn(this, i) result(vgradv)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: i
      complex(kind=dbl)                   :: vgradv(this%jmv)
    end function vgradv_fn

    module subroutine fullnl_sub(this, i)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: i
    end subroutine fullnl_sub
    
    module pure function matica_temp_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),        allocatable  :: matica(:,:)
    end function matica_temp_fn

    module pure function matica_torr_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),        allocatable  :: matica(:,:)
    end function matica_torr_fn

    module pure function matica_mech_fn(this, j_in, a_in) result(matica)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: j_in
      real(kind=dbl),          intent(in) :: a_in
      real(kind=dbl),        allocatable  :: matica(:,:)
    end function matica_mech_fn
    
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
    
    module subroutine solve_temp_sub(this)
      class(T_physicalObject), intent(inout) :: this
    end subroutine solve_temp_sub
  
    module subroutine solve_torr_sub(this)
      class(T_physicalObject), intent(inout) :: this
    end subroutine solve_torr_sub
  
    module subroutine solve_mech_sub(this)
      class(T_physicalObject), intent(inout) :: this
    end subroutine solve_mech_sub
    
    module real(kind=dbl) function nuss_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function nuss_fn

    module real(kind=dbl) function reynolds_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function reynolds_fn

    module real(kind=dbl) function nonzon_reynolds_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function nonzon_reynolds_fn

    module real(kind=dbl) function volume_heating_fn(this)
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

    module real(kind=dbl) function buoyancy_power_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function buoyancy_power_fn

    module real(kind=dbl) function heating_power_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function heating_power_fn

    module real(kind=dbl) function bound_power_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function bound_power_fn

    module real(kind=dbl) function bound_flux_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function bound_flux_fn

    module real(kind=dbl) function advected_heat_fn(this)
      class(T_physicalObject), intent(in) :: this
    end function advected_heat_fn
  end interface
  
end module PhysicalObject