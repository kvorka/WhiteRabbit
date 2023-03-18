module OceanConvMod
  use OceanMod
  implicit none; private
  
  type, extends(T_ocean), public :: T_oceanConv
    contains
    
    procedure, public, pass :: init_sub        => init_oceanConv_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanConv_sub
  
  end type T_oceanConv
  
  contains
  
  subroutine init_oceanConv_sub(this)
    class(T_oceanConv), intent(inout) :: this
    integer                           :: j
    
    call this%init_ocean_sub(); call this%lat_grid%init_vcsv_vcvv_vcvgv_sub()
    
    call this%sol%init_stemp_sub(); call this%sol%init_storr_sub(); call this%sol%init_smech_sub()
    call this%mat%init_mtemp_sub(); call this%mat%init_mtorr_sub(); call this%mat%init_mmech_sub()

    do j=0,this%jmax; call this%mat%temp(j)%fill_sub( matica_temp_fn(this,j,+0.6_dbl), matica_temp_fn(this,j,-0.4_dbl) ); end do
    do j=1,this%jmax; call this%mat%torr(j)%fill_sub( matica_torr_fn(this,j,+0.6_dbl), matica_torr_fn(this,j,-0.4_dbl) ); end do
    do j=1,this%jmax; call this%mat%mech(j)%fill_sub( matica_mech_fn(this,j,+0.6_dbl), matica_mech_fn(this,j,-0.4_dbl) ); end do
    
    allocate( this%nmech(this%jmv,2:this%nd) ) ; this%nmech   = czero
    allocate( this%rmech(this%jmv,2:this%nd) ) ; this%rmech   = czero
    allocate( this%ntemp(this%jms,2:this%nd) ) ; this%ntemp   = czero
    allocate( this%rtemp(this%jms,2:this%nd) ) ; this%rtemp   = czero
    allocate( this%flux_up(this%jms)         ) ; this%flux_up = czero
    
    call this%init_state_sub(); call this%vypis_ocean_sub()
    
  end subroutine init_oceanConv_sub
  
  subroutine time_scheme_oceanConv_sub(this, cf)
    class(T_oceanConv), intent(inout) :: this
    real(kind=dbl),     intent(in)    :: cf
    integer                           :: i, ir1, ir2, ir3, j, jm_int, jml0
    
    !$omp parallel do private (ir1,ir2,ir3,j,jml0) collapse(2)
    do i = 2, this%nd
      do jm_int = 1, this%jms
        j    = this%j_indx(jm_int)
        jml0 = 3*(jm_int-1)
        
        ir1 = 3*(i-1)+1
        ir2 = 6*(i-1)+1
        ir3 = 6*(i-1)+2
        
        if (j == 0) then
          this%rtemp(jm_int,i) = this%mat%temp(j)%multipl_fn(ir1,this%sol%temp(:,jm_int)) + this%ntemp(jm_int,i) * (1-cf)
        else
          this%rtemp(jm_int,i) = this%mat%temp(j)%multipl_fn(ir1,this%sol%temp(:,jm_int)) + this%ntemp(jm_int,i) * (1-cf)
          this%rmech(jml0-1,i) = this%mat%mech(j)%multipl_fn(ir2,this%sol%mech(:,jm_int)) + this%nmech(jml0-1,i) * (1-cf)
          this%rmech(jml0  ,i) = this%mat%torr(j)%multipl_fn(ir1,this%sol%torr(:,jm_int)) + this%nmech(jml0  ,i) * (1-cf)
          this%rmech(jml0+1,i) = this%mat%mech(j)%multipl_fn(ir3,this%sol%mech(:,jm_int)) + this%nmech(jml0+1,i) * (1-cf)
        end if
      end do
    end do
    !$omp end parallel do

    !$omp parallel do
    do i = 2, this%nd
      call fullnl_sub(this, i, this%ntemp(:,i), this%nmech(:,i))
    end do
    !$omp end parallel do
    
    jm_int = 1
      i = 1
        this%sol%temp( 1, 1 ) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
      
      do i = 1, this%nd
        ir1 = 3*(i-1)+1
        
        if (i > 1) then
          this%sol%temp( ir1, 1 ) = this%rtemp(1,i) + cf * this%ntemp(1,i)
        end if

        this%sol%temp( ir1+1 : ir1+2 , 1 ) = czero
      end do
      
      i = this%nd+1
        this%sol%temp( 3*this%nd+1, 1 ) = czero
      
      call this%mat%temp(0)%luSolve_sub( this%sol%temp(:,1) )
      
    !$omp parallel do private (i,ir1,ir2,ir3,j,jml0)
    do jm_int = 2, this%jms
      j    = this%j_indx(jm_int)
      jml0 = 3*(jm_int-1)
      
      i = 1
        this%sol%temp( 1 , jm_int ) = czero
        this%sol%torr( 1 , jm_int ) = czero
        this%sol%mech( 1 , jm_int ) = czero
        this%sol%mech( 2 , jm_int ) = czero
        
      do i = 1, this%nd
        ir1 = 3*(i-1)+1
        ir2 = 6*(i-1)+1
        ir3 = 6*(i-1)+2
        
        if (i > 1) then
          this%sol%temp( ir1 , jm_int ) = this%rtemp( jm_int, i ) + cf * this%ntemp( jm_int, i )
          this%sol%torr( ir1 , jm_int ) = this%rmech( jml0  , i ) + cf * this%nmech( jml0  , i )
          this%sol%mech( ir2 , jm_int ) = this%rmech( jml0-1, i ) + cf * this%nmech( jml0-1, i )
          this%sol%mech( ir3 , jm_int ) = this%rmech( jml0+1, i ) + cf * this%nmech( jml0+1, i )
        end if
        
        this%sol%temp( ir1+1 : ir1+2 , jm_int ) = czero
        this%sol%torr( ir1+1 : ir1+2 , jm_int ) = czero
        this%sol%mech( ir2+2 : ir2+5 , jm_int ) = czero
      end do
      
      i = this%nd+1
        this%sol%temp( 3*this%nd+1 , jm_int ) = czero
        this%sol%torr( 3*this%nd+1 , jm_int ) = czero
        this%sol%mech( 6*this%nd+1 , jm_int ) = czero
        this%sol%mech( 6*this%nd+2 , jm_int ) = czero
        
      call this%mat%temp(j)%luSolve_sub( this%sol%temp(:,jm_int) )
      call this%mat%torr(j)%luSolve_sub( this%sol%torr(:,jm_int) )
      call this%mat%mech(j)%luSolve_sub( this%sol%mech(:,jm_int) )
    end do
    !$omp end parallel do
    
    if (this%mechanic_bnd == 'frees') call this%global_rotation_sub()
    
  end subroutine time_scheme_oceanConv_sub
  
end module OceanConvMod