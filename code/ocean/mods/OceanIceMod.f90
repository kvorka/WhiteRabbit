module OceanIceMod
  use OceanMod
  implicit none

  type, extends(T_ocean), public :: T_oceanice
    real(kind=dbl), private :: ClRoc2

    contains

    procedure, public, pass :: init_sub        => init_oceanice_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanice_sub
  end type T_oceanice

  contains

  subroutine init_oceanice_sub(this)
    class(T_oceanice), intent(inout) :: this
    integer                          :: j
    
    call this%init_ocean_sub()
      call this%lat_grid%init_vcsv_vcvv_vcvgv_sub()

    this%ClRoc2 = Cl_ocean / ( this%Ra * this%Ek**2 / this%Pr )
    
    call this%init_eq_temp_sub()
    call this%init_eq_torr_sub()
    call this%init_eq_mech_sub()
    call this%init_bnd_deformation_sub()

    allocate( this%ntemp(2:this%nd,this%jms) ) ; this%ntemp   = czero
    allocate( this%rtemp(2:this%nd,this%jms) ) ; this%rtemp   = czero
    allocate( this%flux_up(this%jms)         ) ; this%flux_up = czero
    allocate( this%nmech(2:this%nd,this%jmv) ) ; this%nmech   = czero
    allocate( this%rmech(2:this%nd,this%jmv) ) ; this%rmech   = czero
  
    call this%init_state_sub()

  end subroutine init_oceanice_sub

  subroutine time_scheme_oceanice_sub(this, cf)
    class(T_oceanice), intent(inout) :: this
    real(kind=dbl),     intent(in)    :: cf
    integer                           :: i, ir1, ir2, j, jm_int, jml0
    
    !$omp parallel do private (ir1,ir2,j,jml0) collapse(2)
    do jm_int = 1, this%jms
      do i = 2, this%nd
        j    = this%j_indx(jm_int)
        jml0 = 3*(jm_int-1)
        
        ir1 = 3*(i-1)+1
        ir2 = 6*(i-1)+1

        if (j == 0) then
          this%rtemp(i,1) = (1-cf) * this%ntemp(i,1) + this%mat%temp(0)%multipl_fn(ir1,this%sol%temp(:,1))
        else
          this%rtemp(i,jm_int) = (1-cf) * this%ntemp(i,jm_int) + this%mat%temp(j)%multipl_fn(ir1  ,this%sol%temp(:,jm_int))
          this%rmech(i,jml0-1) = (1-cf) * this%nmech(i,jml0-1) + this%mat%mech(j)%multipl_fn(ir2  ,this%sol%mech(:,jm_int))
          this%rmech(i,jml0  ) = (1-cf) * this%nmech(i,jml0  ) + this%mat%torr(j)%multipl_fn(ir1  ,this%sol%torr(:,jm_int))
          this%rmech(i,jml0+1) = (1-cf) * this%nmech(i,jml0+1) + this%mat%mech(j)%multipl_fn(ir2+1,this%sol%mech(:,jm_int))
        end if
      end do
    end do
    !$omp end parallel do

    !$omp parallel do
    do i = 2, this%nd
      call fullnl_sub(this, i, this%ntemp(i,:), this%nmech(i,:))
    end do
    !$omp end parallel do
    
    !$omp parallel workshare
    this%rtemp = this%rtemp + cf * this%ntemp
    this%rmech = this%rmech + cf * this%nmech
    !$omp end parallel workshare
    
    j = 0
      i = 1
        this%sol%temp( 1  , 1 ) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
        this%sol%temp( 2:3, 1 ) = czero
      
      do i = 2, this%nd
        ir1 = 3*(i-1)+1
          this%sol%temp( ir1         , 1 ) = this%rtemp(i,1)
          this%sol%temp( ir1+1:ir1+2 , 1 ) = czero
      end do
      
      i = this%nd+1
        this%sol%temp( 3*this%nd+1, 1 ) = czero
      
      call this%mat%temp(0)%luSolve_sub( this%sol%temp(:,1) ) ; q = real(-this%sol%flux_fn(this%nd,0,0,+1), kind=dbl) / sqrt(4*pi)
      
    !$omp parallel do private (i,ir1,ir2,j,jml0)
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

        if (i > 1) then
          this%sol%temp( ir1   , jm_int ) = this%rtemp( i, jm_int )
          this%sol%torr( ir1   , jm_int ) = this%rmech( i, jml0   )
          this%sol%mech( ir2   , jm_int ) = this%rmech( i, jml0-1 )
          this%sol%mech( ir2+1 , jm_int ) = this%rmech( i, jml0+1 )
        end if
        
        this%sol%temp( ir1+1 : ir1+2 , jm_int ) = czero
        this%sol%torr( ir1+1 : ir1+2 , jm_int ) = czero
        this%sol%mech( ir2+2 : ir2+5 , jm_int ) = czero
      end do
      
      i = this%nd+1
        this%sol%temp( 3*this%nd+1 , jm_int ) = this%ClRoc2 * this%sol%t_up(jm_int) + q * this%sol%u_up(jm_int)
        this%sol%torr( 3*this%nd+1 , jm_int ) = czero
        this%sol%mech( 6*this%nd+1 , jm_int ) = czero
        this%sol%mech( 6*this%nd+2 , jm_int ) = czero
        
      call this%mat%temp(j)%luSolve_sub( this%sol%temp(:,jm_int) )
      call this%mat%torr(j)%luSolve_sub( this%sol%torr(:,jm_int) )
      call this%mat%mech(j)%luSolve_sub( this%sol%mech(:,jm_int) )
    end do
    !$omp end parallel do
    
    if (this%mechanic_bnd == 'frees') call this%global_rotation_sub()

  end subroutine time_scheme_oceanice_sub

end module OceanIceMod