module OceanMod
  use PhysicalObject
  use OceanConstants
  implicit none

  type, extends(T_physicalObject), abstract, public :: T_ocean
    contains
    
    procedure :: init_ocean_sub
    procedure :: speed_sub
    procedure :: set_boundary_deformation_sub
    
    procedure :: vypis_ocean_sub => vypis_ocean_sub
    procedure :: iter_sub        => iter_ocean_sub
    procedure :: init_state_sub  => init_state_ocean_sub
    procedure :: deallocate_sub  => deallocate_ocean_sub

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
    
    call this%init_objects_sub( nd = nd_ocean, jmax = jmax_ocean, r_ud = r_ud_ocean, rgrid = grid_type_ocean, &
                              & gmod = gravity_ocean , g = (1-r_ud_ocean)**2 , noharm = noharm_ocean          )
    
    call this%bnd%init_flux_up_sub()
    
    this%n_iter = n_iter_ocean
    this%cf     = 0.6_dbl
    this%ab     = 1.5_dbl
    
    this%Pr = Pr_ocean
    this%Ra = Ra_ocean
    this%Ek = Ek_ocean
    this%Cl = Kl_ocean
    this%St = St_ocean
    
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
    integer                       :: k, ijm
    real(kind=dbl)                :: avrg_flux
    
    call zero_carray_sub( this%jms, this%bnd%flux_up )
    
    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub()
    end do
    
    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub()
        
        do concurrent ( ijm = 1:this%jms )
          this%bnd%flux_up(ijm) = this%bnd%flux_up(ijm) + this%qr_r_fn(this%nd,ijm)
        end do
    end do
    
    avrg_flux = this%bnd%flux_up(1)%re / s4pi
      do concurrent ( ijm = 1:this%jms )
        this%bnd%flux_up(ijm) = this%bnd%flux_up(ijm) / avrg_flux
      end do
    
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

    write(11,*) this%t, this%dt, this%nuss_fn(), this%reynolds_fn(), this%reynolds_fn(choice='convective')
    write(12,*) this%t, this%dt, this%laws_temp_fn(), this%laws_mech_fn()

    call this%vypis_sub(8, 'data/data_ocean_temp' , 'temperature')
    call this%vypis_sub(8, 'data/data_ocean_veloc', 'velocity'   )
    call this%vypis_sub(8, 'data/data_ocean_flux' , 'flux'       )

    this%poc = this%poc + 1

  end subroutine vypis_ocean_sub
  
  subroutine init_state_ocean_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                           :: i, j, m, ijm, ndI1, jmsI, jmvI
    real(kind=dbl)                    :: ab_help, re, im
    real(kind=dbl),    allocatable    :: r(:)
    complex(kind=dbl), allocatable    :: velc(:), temp(:,:), spher1(:,:), torr(:,:), spher2(:,:)

    if (.not. init_through_file_ocean) then
      do i = 1, this%nd+1
        do j = 0, this%jmax
          do m = 0, j
            ijm = jm(j,m)

            if ((j == 0) .and. (m == 0)) then
              this%sol%temp(3*(i-1)+1,ijm)%re = (this%rad_grid%r(this%nd)/this%rad_grid%rr(i)-1)*this%rad_grid%r(1)*s4pi
            else if (m == 0) then
              call random_number( re )
              this%sol%temp(3*(i-1)+1, ijm)%re = re / 1e3
            else
              call random_number( re ); call random_number( im )
              this%sol%temp(3*(i-1)+1, ijm) = cmplx(re, im, kind=dbl) / 1e3
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

            do ijm = 2, jmsI
              spher1(i,ijm) = velc(3*(ijm-1)-1)
              torr(  i,ijm) = velc(3*(ijm-1)  )
              spher2(i,ijm) = velc(3*(ijm-1)+1)
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
    
    call this%prepare_mat_temp_sub( ijstart=0 , ijend=this%jmax )
    call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
    ab_help = this%ab
    this%ab = one
    
      call this%time_scheme_sub()
      call this%vypis_ocean_sub()
      
    this%ab = ab_help
    
  end subroutine init_state_ocean_sub
  
  subroutine set_boundary_deformation_sub(this, u_up, t_up)
    class(T_ocean),    intent(inout) :: this
    complex(kind=dbl), intent(in)    :: u_up(:), t_up(:)
    integer                          :: jmsmax

    jmsmax = min(size(t_up),this%jms)

    this%bnd%t_up(1:jmsmax) = t_up(1:jmsmax) / this%D_ud
    this%bnd%u_up(1:jmsmax) = u_up(1:jmsmax) / this%D_ud

  end subroutine set_boundary_deformation_sub
  
  subroutine deallocate_ocean_sub(this)
    class(T_ocean), intent(inout) :: this

    close(11); close(12)
    call this%deallocate_objects_sub()

  end subroutine deallocate_ocean_sub
  
end module OceanMod