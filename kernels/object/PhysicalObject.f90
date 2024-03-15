module PhysicalObject
  use Math
  use sph_indexing
  use sph_unitvec_op
  use sph_norms
  use SphericalHarmonics
  use RadialGrid
  use Gravity
  use Solution
  use Matrices
  implicit none
  
  type, abstract, public :: T_physicalObject
    character(len=6)               :: rheology, scaling
    character(len=5)               :: thermal_bnd, mechanic_bnd, grid_type
    logical                        :: noharm, noobj
    integer                        :: nd, jmax, jms, jmv, n_iter, poc
    real(kind=dbl)                 :: t, dt, cf, ab, rd, ru, r_ud, D_ud, gd, gu, Pr, Ra, Ek, Cl, Ds, Raf, Ramu, Rad, Rau
    integer,           allocatable :: j_indx(:)
    complex(kind=dbl), allocatable :: flux_up(:), htide(:,:), rsph1(:,:), rsph2(:,:), rtorr(:,:), rtemp(:,:), &
                                    &                         nsph1(:,:), nsph2(:,:), ntorr(:,:), ntemp(:,:)
    
    type(T_radialGrid)  :: rad_grid
    type(T_lateralGrid) :: lat_grid
    type(T_gravity)     :: gravity
    type(T_matrices)    :: mat
    type(T_solution)    :: sol
    
    contains
    
    procedure, pass :: init_objects_sub       => init_objects_sub
    procedure, pass :: deallocate_objects_sub => deallocate_objects_sub
    
    procedure, pass :: lambda_fn, cp_fn, visc_fn, alpha_fn, set_dt_sub, reynolds_fn, vypis_sub, htide_fn, qr_fn, vr_fn, qr_jm_fn,  &
    & dv_dr_rr_jml_sub, mgradT_rr_jml_sub, buoy_rr_jm_fn, coriolis_vgradv_sub, coriolis_sub, laws_temp_fn, nuss_fn, laws_mech_fn,  &
    & buoy_rr_jml_sub, coriolis_rr_jml_sub, global_rotation_sub, mvgradT_sub, fullnl_sub, matica_temp_fn, matica_mech_fn, vr_jm_fn,&
    & init_eq_temp_sub, init_eq_mech_sub, init_eq_torr_sub, prepare_mat_mech_sub, prepare_mat_temp_sub, prepare_mat_torr_sub,      &
    & solve_temp_sub, solve_torr_sub, solve_mech_sub, volume_heating_fn, laws_force_fn, matica_torr_fn
    
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
    
    module pure complex(kind=dbl) function vr_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function vr_fn
    
    module pure complex(kind=dbl) function qr_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function qr_fn
    
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
    
    module pure complex(kind=dbl) function buoy_rr_jm_fn(this, ir, ijm)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: ir, ijm
    end function buoy_rr_jm_fn
    
    module pure subroutine coriolis_rr_jml_sub(this, v, coriolis)
      class(T_physicalObject), intent(in)    :: this
      complex(kind=dbl),       intent(in)    :: v(:)
      complex(kind=dbl),       intent(inout) :: coriolis(:,:)
    end subroutine coriolis_rr_jml_sub
    
    module pure subroutine buoy_rr_jml_sub(this, ir, T, force)
      class(T_physicalObject), intent(in)    :: this
      integer,                 intent(in)    :: ir
      complex(kind=dbl),       intent(in)    :: T(:)
      complex(kind=dbl),       intent(inout) :: force(:,:)
    end subroutine buoy_rr_jml_sub
    
    module pure subroutine global_rotation_sub(this)
      class(T_physicalObject), intent(inout) :: this
    end subroutine global_rotation_sub
    
    module pure subroutine mvgradT_sub(this, i, mvgradT)
      class(T_physicalObject), intent(in)  :: this
      integer,                 intent(in)  :: i
      complex(kind=dbl),       intent(out) :: mvgradT(:)
    end subroutine mvgradT_sub
    
    module pure subroutine fullnl_sub(this, i)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: i
    end subroutine fullnl_sub

    module pure subroutine coriolis_vgradv_sub(this, i)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: i
    end subroutine coriolis_vgradv_sub

    module pure subroutine coriolis_sub(this, i)
      class(T_physicalObject), intent(inout) :: this
      integer,                 intent(in)    :: i
    end subroutine coriolis_sub
    
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