module OceanIceMod
  use OceanMod
  use MatrixDefinitions
  use NonLinearTerms
  use BalanceEquations
  use VolumeMeassures
  implicit none

  type, extends(T_ocean), public :: T_oceanice
    real(kind=dbl),                 private :: ClRoc2
    complex(kind=dbl), allocatable, private :: nmech(:,:), ntemp(:,:)

    contains

    procedure, public, pass  :: init_sub       => init_oceanice_sub
    procedure, public, pass  :: iter_sub       => iter_oceanice_sub
    procedure, public, pass  :: deallocate_sub => deallocate_oceanice_sub

  end type T_oceanice

  private :: init_oceanice_sub
  private :: init_state_oceanice_sub
  private :: time_scheme_oceanice_sub
  private :: vypis_oceanice_sub
  private :: deallocate_oceanice_sub

  contains

  subroutine init_oceanice_sub(this)
    class(T_oceanice), intent(inout) :: this
    integer                          :: j
    
    call this%init_ocean_sub()

    this%ClRoc2 = Cl_ocean / ( this%Ra * this%Ek**2 / this%Pr )

    call this%lat_grid%init_vcsv_vcvv_vcvgv_sub()
    
    call this%sol%init_layer_u_sub()
    call this%sol%init_stemp_sub(); call this%sol%init_storr_sub(); call this%sol%init_smech_sub()
    call this%mat%init_mtemp_sub(); call this%mat%init_mtorr_sub(); call this%mat%init_mmech_sub()
      do j=0,this%jmax; call this%mat%temp(j)%fill_sub( matica_temp_fn(this,j,+0.6_dbl), matica_temp_fn(this,j,-0.4_dbl)  ); end do
      do j=1,this%jmax; call this%mat%torr(j)%fill_sub( matica_torr_fn(this,j,+0.6_dbl), matica_torr_fn(this,j,-0.4_dbl)  ); end do
      do j=1,this%jmax; call this%mat%mech(j)%fill_sub( matica_mech_fn(this,j,+0.6_dbl), matica_mech_fn(this,j,-0.4_dbl)  ); end do

    allocate( this%nmech(this%jmv,2:this%nd), this%ntemp(this%jms,2:this%nd), this%flux_up(this%jms) )
      this%nmech = cmplx(0._dbl, 0._dbl, kind=dbl); this%ntemp = cmplx(0._dbl, 0._dbl, kind=dbl)
      this%flux_up = cmplx(0._dbl, 0._dbl, kind=dbl)

    open(unit=11, file='data/Nuss.dat', status='new', action='write')
    open(unit=12, file='data/Laws.dat', status='new', action='write')
    
    call init_state_oceanice_sub(this); call vypis_oceanice_sub(this)

  end subroutine iniT_oceanice_sub

  subroutine init_state_oceanice_sub(this)
    class(T_oceanice), intent(inout) :: this
    integer                          :: i, j, m, jm_int, ndI1, jmsI, jmvI
    real(kind=dbl),    allocatable   :: r(:)
    complex(kind=dbl), allocatable   :: velc(:), temp(:,:), spher1(:,:), torr(:,:), spher2(:,:)

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
        spher1 = cmplx(0._dbl, 0._dbl, kind=dbl); spher2 = cmplx(0._dbl, 0._dbl, kind=dbl); torr = cmplx(0._dbl, 0._dbl, kind=dbl)

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
    
    call time_scheme_oceanice_sub(this, cf=1._dbl)

  end subroutine init_state_oceanice_sub

  subroutine iter_oceanice_sub(this, t_bnd, u_bnd)
    class(T_oceanice), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: t_bnd(:), u_bnd(:)
    integer                          :: k

    this%sol%t_up(1:size(t_bnd)) = t_bnd / this%D_ud
    this%sol%u_up(1:size(u_bnd)) = u_bnd / this%D_ud

    this%flux_up = cmplx(0._dbl, 0._dbl, kind=dbl)
    
    do k = 1, 2 * this%n_iter
      call time_scheme_oceanice_sub(this, cf=1.5_dbl)
      if ( k > this%n_iter ) this%flux_up = this%flux_up + this%qr_jm_fn(this%nd)
    end do
    
    this%flux_up = this%flux_up / ( this%flux_up(1)%re / sqrt(4*pi) )

    call vypis_oceanice_sub(this)

  end subroutine iter_oceanice_sub

  subroutine time_scheme_oceanice_sub(this, cf)
    class(T_oceanice), intent(inout) :: this
    real(kind=dbl),    intent(in)    :: cf
    integer                           :: i, jm_int, j, jm1, jm2, jmv1, jmv2
    real(kind=dbl)                   :: q
    complex(kind=dbl)                :: angularMomentum
    complex(kind=dbl), allocatable   :: rtemp(:,:), rmech(:,:)

    this%t = this%t + this%dt

    allocate( rmech(this%jmv,2:this%nd), rtemp(this%jms,2:this%nd) )
      rmech = cmplx(0._dbl, 0._dbl, kind=dbl); rtemp = cmplx(0._dbl, 0._dbl, kind=dbl)

    !$omp parallel do private (j, jm1, jm2, jmv1, jmv2)
      do i = 2, this%nd
        j = 0
          rtemp(1,i) = this%mat%temp(0)%multipl_fn(3*(i-1)+1,this%sol%temp(:,1))
        
        do j = 1, this%jmax
          jm1 = j*(j+1)/2+1 ; jm2 = j*(j+1)/2+j+1 ; jmv1 = 3*(jm1-1)  ; jmv2 = 3*(jm2-1)
        
          rtemp(jm1   :jm2   :1,i) = this%mat%temp(j)%multipl2_fn(3*(i-1)+1,this%sol%temp(:,jm1:jm2))
          rmech(jmv1  :jmv2  :3,i) = this%mat%torr(j)%multipl2_fn(3*(i-1)+1,this%sol%torr(:,jm1:jm2))
          rmech(jmv1-1:jmv2-1:3,i) = this%mat%mech(j)%multipl2_fn(6*(i-1)+1,this%sol%mech(:,jm1:jm2))
          rmech(jmv1+1:jmv2+1:3,i) = this%mat%mech(j)%multipl2_fn(6*(i-1)+2,this%sol%mech(:,jm1:jm2))
        end do
  
        rtemp(:,i) = rtemp(:,i) + this%ntemp(:,i) * (1-cf)
        rmech(:,i) = rmech(:,i) + this%nmech(:,i) * (1-cf)

        call fullnl_sub(this, i, this%ntemp(:,i), this%nmech(:,i))
  
        rtemp(:,i) = rtemp(:,i) + this%ntemp(:,i) * cf
        rmech(:,i) = rmech(:,i) + this%nmech(:,i) * cf
      end do
      !$omp end parallel do

    jm_int = 1
      i = 1
        this%sol%temp( 1   , jm_int ) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
        this%sol%temp( 2:3 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)

      do i = 2, this%nd
        this%sol%temp( 3*(i-1)+1             , jm_int ) = rtemp(jm_int,i)
        this%sol%temp( 3*(i-1)+2 : 3*(i-1)+3 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
      end do

      i = this%nd+1
        this%sol%temp( 3*(i-1)+1 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)

      call this%mat%temp( this%j_indx(jm_int) )%luSolve_sub( this%sol%temp(:,jm_int) )
    
    q = real(-this%sol%flux_fn(this%nd,0,0,+1), kind=dbl) / sqrt(4*pi)

    !$omp parallel do private (i)
    do jm_int = 2, this%jms
      i = 1
        this%sol%temp( 1:3 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
        this%sol%torr( 1:3 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
        this%sol%mech( 1:6 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)

      do i = 2, this%nd
        this%sol%temp( 3*(i-1)+1 , jm_int ) = rtemp(    jm_int     , i )
        this%sol%torr( 3*(i-1)+1 , jm_int ) = rmech( 3*(jm_int-1)  , i )
        this%sol%mech( 6*(i-1)+1 , jm_int ) = rmech( 3*(jm_int-1)-1, i )
        this%sol%mech( 6*(i-1)+2 , jm_int ) = rmech( 3*(jm_int-1)+1, i )
        
        this%sol%temp( 3*(i-1)+2 : 3*(i-1)+3 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
        this%sol%torr( 3*(i-1)+2 : 3*(i-1)+3 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
        this%sol%mech( 6*(i-1)+3 : 6*(i-1)+6 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
      end do

      i = this%nd+1
        this%sol%temp( 3*(i-1)+1             , jm_int ) = this%ClRoc2 * this%sol%t_up(jm_int) + &
                                                        &           q * this%sol%u_up(jm_int)
        this%sol%torr( 3*(i-1)+1             , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
        this%sol%mech( 6*(i-1)+1 : 6*(i-1)+2 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)

      call this%mat%temp( this%j_indx(jm_int) )%luSolve_sub( this%sol%temp(:,jm_int) )
      call this%mat%torr( this%j_indx(jm_int) )%luSolve_sub( this%sol%torr(:,jm_int) )
      call this%mat%mech( this%j_indx(jm_int) )%luSolve_sub( this%sol%mech(:,jm_int) )
    end do
    !$omp end parallel do

    deallocate( rmech, rtemp )

    if (this%mechanic_bnd == 'frees') then
      associate( coeff => ((1/this%r_ud-1)**5) / (1/this%r_ud**5-1) )
        do jm_int = 2, 3
          angularMomentum = 5 * this%rad_grid%intV_fn(this%rad_grid%rr * this%sol%velocity_i_fn(1,jm_int-2,0)) * coeff
          
          do i = 1, this%nd+1
            this%sol%torr(3*(i-1)+1, jm_int) = this%sol%torr(3*(i-1)+1, jm_int) - angularMomentum * this%rad_grid%rr(i)
          end do
        end do
      end associate
    end if

  end subroutine time_scheme_oceanice_sub

  subroutine vypis_oceanice_sub(this)
    !Vypis Nu, Re, zachovania energie a hybnosti, vypis kompletnej diagnostiky oceanu.
    class(T_oceanice), intent(inout) :: this

    write(11,*) this%t, this%dt, nuss_fn(this), reynolds_fn(this), nonzon_reynolds_fn(this)
    write(12,*) this%t, this%dt, laws_temp_fn(this), laws_mech_fn(this)
   
    call this%vypis_sub(8, 'data/data_ocean_temp' , 'temperature')
    call this%vypis_sub(8, 'data/data_ocean_veloc', 'velocity'   )
    call this%vypis_sub(8, 'data/data_ocean_flux' , 'flux'       )
    
    this%poc = this%poc + 1
    
  end subroutine vypis_oceanice_sub

  subroutine deallocate_oceanice_sub(this)
    !Cistenie po vypocte - destruktor pre T_oceanice.
    class(T_oceanice), intent(inout) :: this

    deallocate( this%nmech, this%ntemp, this%flux_up )
    close(11); close(12)
    call this%lat_grid%deallocate_fftw_vcsv_vcvv_vcvgv_sub()
    call this%deallocate_ocean_sub()

  end subroutine deallocate_oceanice_sub

end module OceanIceMod