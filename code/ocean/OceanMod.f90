module OceanMod
  use PhysicalObject
  use OceanConstants
  implicit none

  type, extends(T_physicalObject), abstract, public :: T_ocean
    contains
    
    procedure, pass :: init_ocean_sub
    procedure, pass :: speed_sub
    procedure, pass :: set_boundary_deformation_sub

    procedure, pass :: vypis_ocean_sub => vypis_ocean_sub
    procedure, pass :: iter_sub        => iter_ocean_sub
    procedure, pass :: init_state_sub  => init_state_ocean_sub
    procedure, pass :: deallocate_sub  => deallocate_ocean_sub

    procedure(time_scheme_abstract), deferred, pass :: time_scheme_sub

  end type T_ocean

  abstract interface
    subroutine time_scheme_abstract(this)
       import :: T_ocean, dbl
       class(T_ocean), intent(inout) :: this
    end subroutine time_scheme_abstract
  end interface
  
  contains

  subroutine init_ocean_sub(this)
    class(T_ocean), intent(inout) :: this

    call this%init_objects_sub( nd = nd_ocean, jmax = jmax_ocean, r_ud = r_ud_ocean, &
                              & rgrid = grid_type_ocean, noharm = noharm_ocean       )
    
    call this%gravity%init_sub( gmod = gravity_ocean, g = 1 / this%ru**2 )

    this%n_iter = n_iter_ocean
    this%cf     = 0.6_dbl
    this%ab     = 1.5_dbl

    this%Pr = Pr_ocean
    this%Ra = Ra_ocean
    this%Ek = Ek_ocean
    thiS%Cl = Kl_ocean

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
    
    this%flux_up = czero
    
    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub()
    end do
    
    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub()
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
        call this%time_scheme_sub()
    end do
    
  end subroutine speed_sub

  subroutine vypis_ocean_sub(this)
    class(T_ocean), intent(inout) :: this

    write(11,*) this%t, this%dt, this%nuss_fn(), this%reynolds_fn(), this%nonzon_reynolds_fn()
    write(12,*) this%t, this%dt, this%laws_temp_fn(), this%laws_mech_fn()
     
    call this%vypis_sub(8, 'data/data_ocean_temp' , 'temperature')
    call this%vypis_sub(8, 'data/data_ocean_veloc', 'velocity'   )
    call this%vypis_sub(8, 'data/data_ocean_flux' , 'flux'       )
      
    this%poc = this%poc + 1

  end subroutine vypis_ocean_sub

  subroutine init_state_ocean_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                           :: i, j, m, jm_int, ndI1, jmsI, jmvI
    real(kind=dbl)                    :: ab_help
    real(kind=dbl),    allocatable    :: r(:)
    complex(kind=dbl), allocatable    :: velc(:), temp(:,:), spher1(:,:), torr(:,:), spher2(:,:)

    if (.not. init_through_file_ocean) then
      do i = 1, this%nd+1
        do j = 0, this%jmax
          do m = 0, j
            jm_int = jm(j,m)

            if ((j == 0) .and. (m == 0)) then
              this%sol%temp(3*(i-1)+1,jm_int)%re = (this%rad_grid%r(this%nd)/this%rad_grid%rr(i)-1)*this%rad_grid%r(1)*sqrt(4*pi)
            else if (m == 0) then
              call random_number( this%sol%temp(3*(i-1)+1, jm_int)%re )
              this%sol%temp(3*(i-1)+1, jm_int)%re = this%sol%temp(3*(i-1)+1, jm_int)%re / 1e3
            else
              call random_number( this%sol%temp(3*(i-1)+1, jm_int)%re ); call random_number( this%sol%temp(3*(i-1)+1, jm_int)%im )
              this%sol%temp(3*(i-1)+1, jm_int) = this%sol%temp(3*(i-1)+1, jm_int) / 1e3
            end if

          end do
        end do
      end do

    else
      ndI1 = nd_init_ocean+1; jmsI = jm(jmax_init_ocean,jmax_init_ocean); jmvI = jml(jmax_init_ocean,jmax_init_ocean,+1)

      allocate( r(ndI1), velc(jmvI), temp(ndI1,jmsI), spher1(ndI1,jmsI), spher2(ndI1,jmsI), torr(ndI1,jmsI) )
        spher1 = czero; spher2 = czero; torr = czero

        open(unit=8, file='code/ocean/inittemp', status='old', action='read')
          do i = 1, ndI1
            read(8,*) r(i), temp(i,:)
          end do
        close(8)

        open(unit=8, file='code/ocean/initvelc', status='old', action='read')
          do i = 1, ndI1
            read(8,*) r(i), velc

            do jm_int = 2, jmsI
              spher1(i,jm_int) = velc(3*(jm_int-1)-1)
              torr(  i,jm_int) = velc(3*(jm_int-1)  )
              spher2(i,jm_int) = velc(3*(jm_int-1)+1)
            end do
          end do
        close(8)

      deallocate(velc)

      do i = 1, this%nd+1
        this%sol%temp(3*(i-1)+1,:) = this%rad_grid%interpolation_fn(this%jms, i, r, temp  )
        this%sol%mech(6*(i-1)+1,:) = this%rad_grid%interpolation_fn(this%jms, i, r, spher1)
        this%sol%torr(3*(i-1)+1,:) = this%rad_grid%interpolation_fn(this%jms, i, r, torr  )
        this%sol%mech(6*(i-1)+2,:) = this%rad_grid%interpolation_fn(this%jms, i, r, spher2)
      end do

      deallocate(r, spher1, spher2, torr, temp)
    end if
    
    ab_help = this%ab ; this%ab = 1._dbl
    call this%time_scheme_sub() ; call this%vypis_ocean_sub()
    this%ab = ab_help
    
  end subroutine init_state_ocean_sub

  subroutine set_boundary_deformation_sub(this, u_up, t_up)
    class(T_ocean),    intent(inout) :: this
    complex(kind=dbl), intent(in)    :: u_up(:), t_up(:)
    integer                          :: jmsmax

    jmsmax = min(size(t_up),this%jms)

    this%sol%t_up(1:jmsmax) = t_up(1:jmsmax) / this%D_ud
    this%sol%u_up(1:jmsmax) = u_up(1:jmsmax) / this%D_ud

  end subroutine set_boundary_deformation_sub
  
  subroutine deallocate_ocean_sub(this)
    class(T_ocean), intent(inout) :: this

    close(11); close(12)
    call this%deallocate_objects_sub()

  end subroutine deallocate_ocean_sub

end module OceanMod