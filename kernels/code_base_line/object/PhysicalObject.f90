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
    
    type(T_radialGrid)  :: rad_grid
    type(T_lateralGrid) :: lat_grid
    type(T_gravity)     :: gravity
    type(T_matrices)    :: mat
    type(T_solution)    :: sol
    
    contains
    
    procedure :: init_objects_sub       => init_objects_sub
    procedure :: deallocate_objects_sub => deallocate_objects_sub

    procedure :: lambda_fn
    procedure :: cp_fn
    procedure :: visc_fn
    procedure :: alpha_fn
    procedure :: htide_fn
    procedure :: qr_jm_fn
    procedure :: vr_fn
    procedure :: vr_jm_fn
    procedure :: vypis_sub
    
  end type T_physicalObject
  
  private :: init_objects_sub
  private :: deallocate_objects_sub

  interface
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
  end interface
  
  interface
    module pure complex(kind=dbl) function htide_fn(this, i, jm_int)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: i, jm_int
    end function htide_fn

    module pure complex(kind=dbl) function vr_fn(this, i, j, m)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: i, j, m
    end function vr_fn

    module pure function qr_jm_fn(this, i) result(qr)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: i
      complex(kind=dbl)                   :: qr(this%jms)
    end function qr_jm_fn

    module pure function vr_jm_fn(this, i) result(vr)
      class(T_physicalObject), intent(in) :: this
      integer,                 intent(in) :: i
      complex(kind=dbl)                   :: vr(this%jms)
    end function vr_jm_fn
  end interface

  contains
  
  subroutine init_objects_sub( this, nd, jmax, r_ud, rgrid, noobj, noharm )
    class(T_physicalObject),    intent(inout) :: this
    integer,                    intent(in)    :: nd, jmax
    real(kind=dbl),             intent(in)    :: r_ud
    character(len=*),           intent(in)    :: rgrid
    logical,          optional, intent(in)    :: noobj, noharm
    integer                                   :: j, m
    
    this%nd        = nd
    this%jmax      = jmax
    this%r_ud      = r_ud
    this%rd        = r_ud / ( 1-r_ud )
    this%ru        = 1    / ( 1-r_ud )
    this%grid_type = rgrid

    if (present(noobj)) then
      this%noobj = noobj
    else
      this%noobj = .false.
    end if

    if (present(noharm)) then
      this%noharm = noharm
    else if (this%noobj .eqv. .true.) then
      this%noharm = .true.
    else
      this%noharm = .false.
    end if

    !Pomocne premenne
    this%jms  =   ( this%jmax*(this%jmax+1)/2+this%jmax ) + 1
    this%jmv  = 3*( this%jmax*(this%jmax+1)/2+this%jmax ) + 1

    if (.not. this%noobj) then
      !Inicializuj radialny grid danej velkosti
      call this%rad_grid%init_sub(this%nd, this%r_ud / (1-this%r_ud), 1 / (1-this%r_ud), this%grid_type)
    
      !Inicializuj lateralny grid danej velkosti
      if (.not. this%noharm) call this%lat_grid%init_sub(this%jmax)
    
      !Inicializuj premenne pre riesenie
      call this%sol%init_sub(this%nd, this%jmax)
    
      !Inicializuj premenne pre matice
      call this%mat%init_sub(this%nd, this%jmax, this%grid_type)

      !Indexacia
      allocate( this%j_indx(this%jms) )
      do j = 0, this%jmax
        do m = 0, j
          this%j_indx(j*(j+1)/2+m+1) = j
        end do
      end do
    end if
      
    !Cas a casovy krok, vypis
    this%poc = 0
    this%t   = 0._dbl

    if (this%noobj) then
      call this%rad_grid%init_sub(this%nd, this%r_ud / (1-this%r_ud), 1 / (1-this%r_ud), this%grid_type)
      this%dt  = 0.49_dbl * ( this%rad_grid%r(2)-this%rad_grid%r(1) )**2
      call this%rad_grid%deallocate_sub()
    else
      this%dt  = 0.49_dbl * ( this%rad_grid%r(2)-this%rad_grid%r(1) )**2
    end if
    
  end subroutine init_objects_sub
  
  subroutine vypis_sub(this, filenum, path, quantity)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: filenum
    character(len=*),        intent(in) :: path, quantity
    integer                             :: i, j, m

    select case (quantity)
      case ('temperature')
        open(unit=filenum, file=path//'/Temp-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do i = 1, this%nd+1
            write(filenum,*) this%rad_grid%rr(i), this%sol%temp_jm_fn(i)
          end do
        close(filenum)

      case ('velocity')
        open(unit=filenum, file=path//'/Velc-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do i = 1, this%nd+1
            write(filenum,*) this%rad_grid%rr(i), this%sol%velocity_jml_fn(i)
          end do
        close(filenum)
      
      case ('flux')
        open(unit=filenum, file=path//'/Flux-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do j = 0, this%jmax
            do m = 0, j
              write(filenum,*) j, m, this%flux_up(jm(j,m))
            end do
          end do
        close(filenum)
      
      case ('topo')
        open(unit=filenum, file=path//'/Topo-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do j = 1, this%jmax
            do m = 0, j
              write(filenum,*) j, m, this%sol%t_dn(jm(j,m)) * this%D_ud, this%sol%t_up(jm(j,m)) * this%D_ud
            end do
          end do
        close(filenum)
      
      case('shape')
        open(unit=7, file=path//'/Shape-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do j = 1, this%jmax
            do m = 0, j
              write(7,*) j, m, this%sol%u_dn(jm(j,m)) * this%D_ud, this%sol%u_up(jm(j,m)) * this%D_ud
            end do
          end do
        close(7)
    end select

  end subroutine vypis_sub

  subroutine deallocate_objects_sub(this)
    class(T_physicalObject), intent(inout) :: this
    
    if (.not. this%noobj) then
      deallocate( this%j_indx )
      
      if (.not. this%noharm) call this%lat_grid%deallocate_sub()
      call this%rad_grid%deallocate_sub()
      call this%mat%deallocate_sub()
      call this%sol%deallocate_sub()
    end if
    
  end subroutine deallocate_objects_sub
  
end module PhysicalObject