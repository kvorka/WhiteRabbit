module OceanConvMod
  use OceanMod
  implicit none
  private
  
  type, extends(T_ocean), public :: T_oceanConv
    complex(kind=dbl), allocatable, private :: nmech(:,:), ntemp(:,:)
    
    contains
    
    procedure, public, pass :: init_sub        => init_oceanConv_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanConv_sub
    procedure, public, pass :: deallocate_sub  => deallocate_oceanConv_sub
  
  end type T_oceanConv
  
  contains
  
  subroutine init_oceanConv_sub(this)
    class(T_oceanConv), intent(inout) :: this
    integer                           :: j
    
    call this%init_ocean_sub()
    call this%lat_grid%init_vcsv_vcvv_vcvgv_sub()
    
    call this%sol%init_stemp_sub(); call this%sol%init_storr_sub(); call this%sol%init_smech_sub()
    call this%mat%init_mtemp_sub(); call this%mat%init_mtorr_sub(); call this%mat%init_mmech_sub()

    do j=0,this%jmax; call this%mat%temp(j)%fill_sub( matica_temp_fn(this,j,+0.6_dbl), matica_temp_fn(this,j,-0.4_dbl)  ); end do
    do j=1,this%jmax; call this%mat%torr(j)%fill_sub( matica_torr_fn(this,j,+0.6_dbl), matica_torr_fn(this,j,-0.4_dbl)  ); end do
    do j=1,this%jmax; call this%mat%mech(j)%fill_sub( matica_mech_fn(this,j,+0.6_dbl), matica_mech_fn(this,j,-0.4_dbl)  ); end do
    
    allocate( this%nmech(this%jmv,2:this%nd), this%ntemp(this%jms,2:this%nd), this%flux_up(this%jms) )
      this%nmech = cmplx(0._dbl, 0._dbl, kind=dbl); this%ntemp = cmplx(0._dbl, 0._dbl, kind=dbl)
      this%flux_up = cmplx(0._dbl, 0._dbl, kind=dbl)
    
    call this%init_state_sub(); call this%vypis_ocean_sub()
    
  end subroutine init_oceanConv_sub
  
  subroutine time_scheme_oceanConv_sub(this, cf)
    class(T_oceanConv), intent(inout) :: this
    real(kind=dbl),     intent(in)    :: cf
    integer                           :: i, jm_int, j, jm1, jm2, jmv1, jmv2
    complex(kind=dbl)                 :: angularMomentum
    complex(kind=dbl), allocatable    :: rtemp(:,:), rmech(:,:)
    
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
          this%sol%temp( 3*(i-1)+1             , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
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
    
  end subroutine time_scheme_oceanConv_sub
  
  subroutine deallocate_oceanConv_sub(this)
    class(T_oceanConv), intent(inout) :: this
    
    close(11); close(12)
    deallocate( this%nmech, this%ntemp )

    call this%lat_grid%deallocate_fftw_vcsv_vcvv_vcvgv_sub()
    call this%deallocate_ocean_sub()
  
  end subroutine deallocate_oceanConv_sub
  
end module OceanConvMod
