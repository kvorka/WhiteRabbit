module OceanTidesMod
  use OceanMod
  use MatrixDefinitions
  use NonLinearTerms
  use BalanceEquations
  use VolumeMeassures
  implicit none

  type, extends(T_ocean), public :: T_oceanTides
    complex(kind=dbl), allocatable, private :: nmech(:,:), v201(:), v203(:), v221(:), v223(:)
    real(kind=dbl)                          :: heating
    integer                                 :: number_of_periods

    contains

    procedure, public, pass  :: init_sub       => init_oceanTides_sub
    procedure, public, pass  :: iter_sub       => iter_oceanTides_sub
    procedure, public, pass  :: deallocate_sub => deallocate_oceanTides_sub

  end type T_oceanTides

  private :: init_oceanTides_sub
  private :: iter_oceanTides_sub
  private :: time_scheme_oceanTides_sub
  private :: vypis_oceanTides_sub
  private :: deallocate_oceanTides_sub

  contains

  subroutine init_oceanTides_sub(this)
    class(T_oceanTides),  intent(inout) :: this
    integer                             :: j, i
    complex(kind=dbl), allocatable      :: u201(:), u203(:), u221(:), u223(:)
    
    call this%init_ocean_sub()
    
    this%dt = 2 * pi / this%n_iter
    this%number_of_periods = 0

    call this%sol%init_storr_sub(); call this%sol%init_smech_sub()
    call this%mat%init_mtorr_sub(); call this%mat%init_mmech_sub()
      do j=1,this%jmax; call this%mat%torr(j)%fill_sub( matica_torr_fn(this,j,+1._dbl), matica_torr_fn(this,j,0._dbl)  ); end do 
      do j=1,this%jmax; call this%mat%mech(j)%fill_sub( matica_mech_fn(this,j,+1._dbl), matica_mech_fn(this,j,0._dbl)  ); end do

    allocate( this%nmech(this%jmv,2:this%nd) )
      this%nmech = cmplx(0._dbl, 0._dbl, kind=dbl)

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

    open(unit=1, file='data/output_ocean_heating', status='old', action='write')
    
  end subroutine init_oceanTides_sub

  subroutine iter_oceanTides_sub(this)
    class(T_oceanTides), intent(inout) :: this
    integer                            :: k
    
    this%heating = 0._dbl; this%number_of_periods = this%number_of_periods + 1

    do k = 1, this%n_iter
      call time_scheme_oceanTides_sub(this, cf=1.5_dbl, k_per=k)
      this%heating = this%heating + volume_heating_fn(this)
    end do

    this%heating = this%heating / this%n_iter

    call vypis_oceanTides_sub(this)

  end subroutine iter_oceanTides_sub

  subroutine time_scheme_oceanTides_sub(this, cf, k_per)
    class(T_oceanTides),  intent(inout) :: this
    real(kind=dbl),       intent(in)    :: cf
    integer,              intent(in)    :: k_per
    integer                             :: i, jm_int
    complex(kind=dbl),    allocatable   :: rmech(:,:)

    allocate( rmech(this%jmv,2:this%nd) ); rmech = cmplx(0._dbl, 0._dbl, kind=dbl)

    !$omp parallel do private (jm_int)
    do i = 2, this%nd
      do jm_int = 2, this%jms
        rmech(3*(jm_int-1)-1,i) = this%mat%mech( this%j_indx(jm_int) )%multipl_fn( 6*(i-1)+1, this%sol%mech(:,jm_int) )
        rmech(3*(jm_int-1)  ,i) = this%mat%torr( this%j_indx(jm_int) )%multipl_fn( 3*(i-1)+1, this%sol%torr(:,jm_int) )
        rmech(3*(jm_int-1)+1,i) = this%mat%mech( this%j_indx(jm_int) )%multipl_fn( 6*(i-1)+2, this%sol%mech(:,jm_int) )
      end do

      rmech(:,i) = rmech(:,i) + this%nmech(:,i) * (1-cf)

      if (.not. this%noharm) then
        this%nmech(:,i) = 2 * coriolis_fn(this,i) + vgradv_fn(this,i)
      else
        this%nmech(:,i) = 2 * coriolis_fn(this,i)
      end if

      rmech(:,i) = rmech(:,i) + this%nmech(:,i) * (  cf)
    end do
    !$omp end parallel do

    !$omp parallel do private (i)
    do jm_int = 2, this%jms
      i = 1
        this%sol%torr( 1:3 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
        this%sol%mech( 1:6 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)

      do i = 2, this%nd
        this%sol%torr( 3*(i-1)+1 , jm_int ) = rmech( 3*(jm_int-1)  , i )
        this%sol%mech( 6*(i-1)+1 , jm_int ) = rmech( 3*(jm_int-1)-1, i )
        this%sol%mech( 6*(i-1)+2 , jm_int ) = rmech( 3*(jm_int-1)+1, i )
        
        this%sol%torr( 3*(i-1)+2 : 3*(i-1)+3 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
        this%sol%mech( 6*(i-1)+3 : 6*(i-1)+6 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
      end do

      i = this%nd+1
        this%sol%torr( 3*(i-1)+1 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)

        if (jm_int == 4) then
          this%sol%mech( 6*(i-1)+1 , jm_int ) = this%v201(k_per)
          this%sol%mech( 6*(i-1)+2 , jm_int ) = this%v203(k_per)
        else if (jm_int == 6) then
          this%sol%mech( 6*(i-1)+1 , jm_int ) = this%v221(k_per)
          this%sol%mech( 6*(i-1)+2 , jm_int ) = this%v223(k_per)
        else
          this%sol%mech( 6*(i-1)+1 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
          this%sol%mech( 6*(i-1)+2 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
        end if
        
      call this%mat%torr( this%j_indx(jm_int) )%luSolve_sub( this%sol%torr(:,jm_int) )
      call this%mat%mech( this%j_indx(jm_int) )%luSolve_sub( this%sol%mech(:,jm_int) )
    end do
    !$omp end parallel do

    deallocate( rmech )

  end subroutine time_scheme_oceanTides_sub

  subroutine vypis_oceanTides_sub(this)
    class(T_oceanTides), intent(inout) :: this

    write(1,*) this%number_of_periods, stress_dim * this%heating / 1e6
    write(*,*) this%number_of_periods, stress_dim * this%heating / 1e6
    
  end subroutine vypis_oceanTides_sub

  subroutine deallocate_oceanTides_sub(this)
    class(T_oceanTides), intent(inout) :: this

    close(1); deallocate( this%nmech )
    call this%deallocate_ocean_sub()

  end subroutine deallocate_oceanTides_sub

end module OceanTidesMod