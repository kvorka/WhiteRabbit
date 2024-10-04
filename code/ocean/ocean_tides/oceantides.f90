module oceantides
  use ocean
  implicit none

  type, extends(T_ocean), public :: T_oceanTides
    complex(kind=dbl), allocatable, private :: v201(:), v203(:), v221(:), v223(:)
    real(kind=dbl)                          :: heating
    integer                                 :: number_of_periods, k_of_period

    contains

    procedure, public, pass :: init_sub        => init_oceanTides_sub
    procedure, public, pass :: iter_sub        => iter_oceanTides_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanTides_sub
    procedure, public, pass :: vypis_ocean_sub => vypis_oceanTides_sub

  end type T_oceanTides

  contains

  subroutine init_oceanTides_sub(this)
    class(T_oceanTides),  intent(inout) :: this
    integer                             :: i
    complex(kind=dbl), allocatable      :: u201(:), u203(:), u221(:), u223(:)
    
    call this%init_ocean_sub()
    
    this%dt = 2 * pi / this%n_iter ; this%number_of_periods = 0
    
    call this%init_eq_torr_sub()
      allocate( this%ntorr(this%jms,2:this%nd) )
      this%ntorr = czero
      
      call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    
    call this%init_eq_mech_sub()
      allocate( this%nsph1(this%jms,2:this%nd), this%nsph2(this%jms,2:this%nd) )
      this%nsph1 = czero
      this%nsph2 = czero
      
      call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
    allocate( this%v201(this%n_iter), this%v203(this%n_iter), this%v221(this%n_iter), this%v223(this%n_iter) )
      allocate( u201(this%n_iter), u203(this%n_iter), u221(this%n_iter), u223(this%n_iter) )
        open(unit=1, file='code/ocean/OceanTides/upper_bound_disp.txt', status='old', action='read')
          do
            read(1,*) i, u201(i), u203(i), u221(i), u223(i)
            if (i == this%n_iter) exit
          end do
        close(1)

        u201(:) = u201(:) / D_ud_ocean; u203(:) = u203(:) / D_ud_ocean
        u221(:) = u221(:) / D_ud_ocean; u223(:) = u223(:) / D_ud_ocean

        this%v201(1) = ( u201(1) - u201(this%n_iter) ) / this%dt; this%v203(1) = ( u203(1) - u203(this%n_iter) ) / this%dt
        this%v221(1) = ( u221(1) - u221(this%n_iter) ) / this%dt; this%v223(1) = ( u223(1) - u223(this%n_iter) ) / this%dt
        
        do i = 2, this%n_iter
          this%v201(i) = ( u201(i) - u201(i-1) ) / this%dt; this%v203(i) = ( u203(i) - u203(i-1) ) / this%dt
          this%v221(i) = ( u221(i) - u221(i-1) ) / this%dt; this%v223(i) = ( u223(i) - u223(i-1) ) / this%dt
        end do

      deallocate( u201, u203, u221, u223 )
    
  end subroutine init_oceanTides_sub

  subroutine iter_oceanTides_sub(this)
    class(T_oceanTides), intent(inout) :: this
    integer                            :: k
    
    this%heating = zero; this%number_of_periods = this%number_of_periods + 1

    do k = 1, this%n_iter
      this%k_of_period = k ; call this%time_scheme_sub()
      this%heating = this%heating + this%viscdissip_power_fn()
    end do

    this%heating = this%heating / this%n_iter

    call vypis_oceanTides_sub(this)

  end subroutine iter_oceanTides_sub

  subroutine time_scheme_oceanTides_sub(this)
    class(T_oceanTides),  intent(inout) :: this
    integer                             :: ir, ijm
    
    !$omp parallel
    !$omp do collapse (2)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtorr(ir,ijm) = (1-this%ab) * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = (1-this%ab) * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = (1-this%ab) * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end do
    
    if ( this%noharm ) then
      !$omp parallel do
      do ir = 2, this%nd
        call this%coriolis_sub(ir)
      end do
      !$omp end parallel do
    else
      !$omp parallel do
      do ir = 2, this%nd
        call this%coriolis_vgradv_sub(ir)
      end do
      !$omp end parallel do
    end if
    
    !$omp do collapse (2)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtorr(ir,ijm) = this%rtorr(ir,ijm) + this%ab * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + this%ab * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + this%ab * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end do
    !$omp end parallel
    
    !$omp parallel do
    do ijm = 2, this%jms
      this%rtorr(1,ijm) = czero
      this%rsph1(1,ijm) = czero
      this%rsph2(1,ijm) = czero
      
      this%rtorr(this%nd+1,ijm) = czero
      
      if (ijm == 4) then
        this%rsph1(this%nd+1,ijm) = this%v201(this%k_of_period)
        this%rsph2(this%nd+1,ijm) = this%v203(this%k_of_period)
      else if (ijm == 6) then
        this%rsph1(this%nd+1,ijm) = this%v221(this%k_of_period)
        this%rsph2(this%nd+1,ijm) = this%v223(this%k_of_period)
      else
        this%rsph1(this%nd+1,ijm) = czero
        this%rsph2(this%nd+1,ijm) = czero
      end if
    end do
    !$omp end parallel do
    
    call this%solve_torr_sub( ijmstart=2 , ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.true. )
    call this%solve_mech_sub( ijmstart=2 , ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.true. )
    
  end subroutine time_scheme_oceanTides_sub
  
  subroutine vypis_oceanTides_sub(this)
    class(T_oceanTides), intent(inout) :: this

    write(11,*) this%number_of_periods, stress_dim * this%heating / 1e6
    
  end subroutine vypis_oceanTides_sub
  
end module oceantides