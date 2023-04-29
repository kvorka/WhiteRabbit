module OceanMod
  use Math
  use PhysicalObject
  use OceanConstants
  use VolumeMeassures
  use BalanceEquations
  use MatrixDefinitions
  use NonLinearTerms
  implicit none

  type, extends(T_physicalObject), abstract, public :: T_ocean
    complex(kind=dbl), allocatable :: rsph1(:,:), rsph2(:,:), rtorr(:,:), rtemp(:,:) 
    complex(kind=dbl), allocatable :: nsph1(:,:), nsph2(:,:), ntorr(:,:), ntemp(:,:)

    contains
    
    procedure, pass :: init_ocean_sub
    procedure, pass :: speed_sub

    procedure, pass :: set_boundary_deformation_sub
    procedure, pass :: global_rotation_sub

    procedure, pass :: init_eq_temp_sub
    procedure, pass :: init_eq_mech_sub
    procedure, pass :: init_eq_torr_sub
    procedure, pass :: init_bnd_deformation_sub

    procedure, pass :: solve_temp_sub
    procedure, pass :: solve_torr_sub
    procedure, pass :: solve_mech_sub

    procedure, pass :: vypis_ocean_sub => vypis_ocean_sub
    procedure, pass :: iter_sub        => iter_ocean_sub
    procedure, pass :: init_state_sub  => init_state_ocean_sub
    procedure, pass :: deallocate_sub  => deallocate_ocean_sub

    procedure(time_scheme_abstract), deferred, pass :: time_scheme_sub

  end type T_ocean

  abstract interface
    subroutine time_scheme_abstract(this, cf)
       import :: T_ocean, dbl
       class(T_ocean), intent(inout) :: this
       real(kind=dbl), intent(in)    :: cf
    end subroutine time_scheme_abstract
  end interface
  
  interface
    module subroutine init_eq_temp_sub(this,rhs,nl)
      class(T_ocean), intent(inout) :: this
      logical,        intent(in)    :: rhs, nl
    end subroutine init_eq_temp_sub

    module subroutine init_eq_torr_sub(this,rhs,nl)
      class(T_ocean), intent(inout) :: this
      logical,        intent(in)    :: rhs, nl
    end subroutine init_eq_torr_sub

    module subroutine init_eq_mech_sub(this,rhs,nl)
      class(T_ocean), intent(inout) :: this
      logical,        intent(in)    :: rhs, nl
    end subroutine init_eq_mech_sub

    module subroutine init_bnd_deformation_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine init_bnd_deformation_sub

    module subroutine set_boundary_deformation_sub(this, u_up, t_up)
      class(T_ocean),    intent(inout) :: this
      complex(kind=dbl), intent(in)    :: u_up(:), t_up(:)
    end subroutine set_boundary_deformation_sub

    module subroutine global_rotation_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine global_rotation_sub

    module subroutine init_state_ocean_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine init_state_ocean_sub
  end interface

  interface
    module subroutine solve_temp_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine solve_temp_sub

    module subroutine solve_torr_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine solve_torr_sub

    module subroutine solve_mech_sub(this)
      class(T_ocean), intent(inout) :: this
    end subroutine solve_mech_sub
  end interface
  
  contains

  subroutine init_ocean_sub(this)
    class(T_ocean), intent(inout) :: this

    call this%init_objects_sub( nd = nd_ocean, jmax = jmax_ocean, r_ud = r_ud_ocean, &
                              & rgrid = grid_type_ocean, noharm = noharm_ocean       )
    
    call this%gravity%init_sub( gmod = gravity_ocean, g = 1 / this%ru**2 )

    this%n_iter = n_iter_ocean

    this%Pr = Pr_ocean
    this%Ra = Ra_ocean
    this%Ek = Ek_ocean

    this%D_ud         = D_ud_ocean
    this%rheology     = rheology_ocean
    this%mechanic_bnd = mechanic_bnd_ocean
    this%thermal_bnd  = thermal_bnd_ocean
    this%scaling      = scaling_ocean

    open(unit=11, file='data/Nuss.dat', status='new', action='write')
    open(unit=12, file='data/Laws.dat', status='new', action='write')

  end subroutine init_ocean_sub

  subroutine iter_ocean_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: k

    this%flux_up = cmplx(0._dbl, 0._dbl, kind=dbl)
    
    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub(cf=1.5_dbl)
    end do

    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub(cf=1.5_dbl)
        this%flux_up = this%flux_up + this%qr_jm_fn(this%nd)
    end do
    
    this%flux_up = this%flux_up / ( this%flux_up(1)%re / sqrt(4*pi) )

    call this%vypis_ocean_sub()

  end subroutine iter_ocean_sub

  subroutine speed_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: k
    
    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub(cf=1.5_dbl)
    end do
    
  end subroutine speed_sub

  subroutine vypis_ocean_sub(this)
    class(T_ocean), intent(inout) :: this

    write(11,*) this%t, this%dt, nuss_fn(this), reynolds_fn(this), nonzon_reynolds_fn(this)
    write(12,*) this%t, this%dt, laws_temp_fn(this), laws_mech_fn(this)
     
    call this%vypis_sub(8, 'data/data_ocean_temp' , 'temperature')
    call this%vypis_sub(8, 'data/data_ocean_veloc', 'velocity'   )
    call this%vypis_sub(8, 'data/data_ocean_flux' , 'flux'       )
      
    this%poc = this%poc + 1

  end subroutine vypis_ocean_sub

  subroutine deallocate_ocean_sub(this)
    class(T_ocean), intent(inout) :: this

    close(11); close(12)
    
    if ( allocated(this%flux_up) ) deallocate( this%flux_up )

    if ( allocated(this%ntemp) ) deallocate( this%ntemp )
    if ( allocated(this%ntorr) ) deallocate( this%ntorr )
    if ( allocated(this%nsph1) ) deallocate( this%nsph1 )
    if ( allocated(this%nsph2) ) deallocate( this%nsph2 )

    if ( allocated(this%rtemp) ) deallocate( this%rtemp )
    if ( allocated(this%rtorr) ) deallocate( this%rtorr )
    if ( allocated(this%rsph1) ) deallocate( this%rsph1 )
    if ( allocated(this%rsph2) ) deallocate( this%rsph2 )

    call this%deallocate_objects_sub()

  end subroutine deallocate_ocean_sub

end module OceanMod