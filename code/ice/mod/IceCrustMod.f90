module IceCrustMod
  use IceMod
  use IceTidesMod
  implicit none
  
  type, extends(T_ice), public :: T_iceCrust
    contains
    
    procedure, public, pass :: init_sub        => init_iceCrust_sub
    procedure, public, pass :: time_scheme_sub => iter_iceCrust_sub
    
    procedure, private, pass :: EE_sub        => EE_iceCrust_sub
    procedure, private, pass :: EE_temp00_sub => EE_temp00_iceCrust_sub
    procedure, private, pass :: EE_temp_sub   => EE_temp_iceCrust_sub
    procedure, private, pass :: EE_mech_sub   => EE_mech_iceCrust_sub
    procedure, private, pass :: EE_topo_sub   => EE_topo_iceCrust_sub
    
  end type T_iceCrust
  
  private :: init_iceCrust_sub
  private :: vypis_iceCrust_sub
  
  private :: EE_iceCrust_sub
  private :: EE_temp00_iceCrust_sub
  private :: EE_temp_iceCrust_sub
  private :: EE_mech_iceCrust_sub
  private :: EE_topo_iceCrust_sub
  
  contains
  
  subroutine init_iceCrust_sub(this, notides)
    class(T_iceCrust), intent(inout) :: this
    logical, optional, intent(in)    :: notides
    type(T_iceTides)                 :: tides
    integer                          :: i, j, m
    real(kind=dbl)                   :: radius
    complex(kind=dbl), allocatable   :: flux_bnd(:)
   
    call this%init_ice_sub(jmax_in = jmax_ice, rheol_in = 'viscos', n_iter = n_iter_ice)
    call this%lat_grid%init_vcvv_sub()
    
    this%cf = 1._dbl

    call this%init_eq_temp_sub( rhs=.true. , nl=.true.  )
    call this%init_eq_mech_sub( rhs=.true. , nl=.false. )
    
    call tides%init_sub()
    call tides%deallocate_sub()
    write(*,*) 'OK'
    
    !if ( notides ) then
    !  open(unit=1, file='data/data_ice_tides/Temp_tides.dat', status='old', action='read')
    !    do i = 1, this%nd+1; read(1,*) radius, this%sol%temp(3*(i-1)+1,1); end do
    !  close(1)
    !  
    !else
      open(unit=1, file='data/data_ice_tides/Temp_tides.dat', status='old', action='read')
        do i = 1, this%nd+1; read(1,*) radius, this%sol%temp(3*(i-1)+1,1); end do
      close(1)
    
      !open(unit=1, file='data/data_ice_tides/tidal_heating.dat', status='old', action='read')
      !  do i = 1, this%nd; read(1,*) radius, this%htide(i,:); end do
      !close(1)
    !end if
    
    allocate( flux_bnd(this%jms) ); flux_bnd = czero
      
    do
      call this%EE_sub(flux_bnd)
      !write(*,*) c2r_fn(this%sol%u_up(4)) * this%D_ud , c2r_fn(this%htide_fn(3,4)), c2r_fn(this%ntemp(4,3)), &
      !          & abs(this%sol%v_up(4) * this%dt / this%sol%u_up(4))
      if ( abs(this%sol%v_up(4) * this%dt / this%sol%u_up(4)) < 1e-9 ) exit
    end do

    write(*,*) c2r_fn(this%qr_fn(1,4)) / ( c2r_fn(this%qr_fn(1,1)) / sqrt(4*pi) )
    
    deallocate( flux_bnd )
    
    write(*,*) 'Icy crust quasi hydrostatic' ; call vypis_iceCrust_sub(this)

  end subroutine init_iceCrust_sub

  subroutine vypis_iceCrust_sub(this)
    class(T_iceCrust), intent(inout) :: this
      
    call this%vypis_sub(8, 'data/data_ice_topo' , 'topo' )
    call this%vypis_sub(8, 'data/data_ice_shape', 'shape')
    
    this%poc = this%poc + 1
    
  end subroutine vypis_iceCrust_sub
  
  subroutine iter_iceCrust_sub(this, flux_bnd)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux_bnd(:)
    integer                          :: n
    complex(kind=dbl), allocatable   :: flux(:)
    
    flux = flux_bnd(1:this%jms)
    
    do n = 1, this%n_iter
      call this%EE_sub( flux )
    end do

    call vypis_iceCrust_sub(this)
    
    deallocate( flux )
    
  end subroutine iter_iceCrust_sub
  
  subroutine EE_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(inout) :: flux(:)
    integer                          :: ir
    
    this%t = this%t + this%dt
    
    !$omp parallel do
    do ir = 2, this%nd
      this%ntemp(:,ir) = -this%vgradT_fn(ir)
    end do
    !$omp end parallel do
    
    call this%EE_temp00_sub(flux)
    call this%EE_temp_sub(flux)
    call this%EE_mech_sub(flux)
    call this%EE_topo_sub(flux)
    
  end subroutine EE_iceCrust_sub

  subroutine EE_temp00_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(inout) :: flux(:)
    integer                          :: ir, is
    complex(kind=dbl), allocatable   :: Temp(:), Temp1(:)

    allocate( Temp(this%nd+1), Temp1(this%nd+1) ); Temp = this%sol%temp_i_fn(1)
    
    do
      Temp1 = this%sol%temp_i_fn(1)
      call this%mat%temp(0)%fill_sub( this%matica_temp_fn(j_in=0, a_in=  this%cf), &
                                    & this%matica_temp_fn(j_in=0, a_in=1-this%cf)  )
    
      ir = 1
        this%sol%temp(1,1) = cs4pi
        this%sol%temp(2,1) = czero
        this%sol%temp(3,1) = czero
  
      do ir = 2, this%nd
        is = 3*(ir-1)+1
        
        this%sol%temp(is  ,1) = Temp(ir) / this%dt + this%htide_fn(ir,1) + this%ntemp(1,ir)
        this%sol%temp(is+1,1) = czero
        this%sol%temp(is+2,1) = czero
      end do
    
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,1) = czero
    
      call this%mat%temp(0)%luSolve_sub(this%sol%temp(:,1))
    
      if ( maxval(abs(this%sol%temp_i_fn(1) - Temp1)/abs(Temp1)) < 1e-5 ) exit       
    end do
    
    deallocate( Temp, Temp1 )
    
    flux = flux * c2r_fn( -this%sol%flux_fn(1,1,1) / sqrt(4*pi) )
    
  end subroutine EE_temp00_iceCrust_sub
  
  subroutine EE_temp_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux(:)
    integer                          :: ir, ijm
    
    !$omp parallel do private (ir)
    do ijm = 2, this%jms
      ir = 1
        this%rtemp(1,ijm) = -( this%sol%u_dn(ijm) + ( this%vr_fn(1,ijm) + this%Raf * flux(ijm) ) * this%dt )
      
      do concurrent ( ir = 2:this%nd )
        this%rtemp(ir,ijm) = this%htide_fn(ir,ijm) + this%ntemp(ijm,ir)
      end do
      
      ir = this%nd+1
        this%rtemp(this%nd+1,ijm) = -( this%sol%u_up(ijm) + this%vr_fn(this%nd,ijm) * this%dt )
    end do
    !$omp end parallel do
    
    call this%solve_temp_sub( ijmstart=2, rematrix=.true. )
    
  end subroutine EE_temp_iceCrust_sub
  
  subroutine EE_mech_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux(:)
    integer                          :: ir, ij, ijm
    
    !$omp parallel do private (ir,ij)
    do ijm = 2, this%jms
      ij = this%j_indx(ijm)

      ir = 1
        this%rsph1(ir,ijm) = -( this%sol%u_dn(ijm) - this%Vdelta_fn(1,ijm) - this%Raf * (this%qr_fn(ir,ijm) - flux(ijm)) * this%dt )
        this%rsph2(ir,ijm) = czero
      
      do concurrent ( ir = 2:this%nd )
        this%rsph1(ir,ijm) = -sqrt((ij  )/(2*ij+1._dbl)) * this%buoy_rr_fn(ir,ijm)
        this%rsph2(ir,ijm) = +sqrt((ij+1)/(2*ij+1._dbl)) * this%buoy_rr_fn(ir,ijm)
      end do
      
      ir = this%nd+1
        this%rsph1(ir,ijm) = czero
        this%rsph2(ir,ijm) = -( this%sol%u_up(ijm) - this%Vdelta_fn(this%nd,ijm) )
    end do
    !$omp end parallel do

    call this%solve_mech_sub( ijmstart=2, rematrix=.true. )
    
  end subroutine EE_mech_iceCrust_sub
  
  subroutine EE_topo_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux(:)
    integer                          :: ijm
    
    this%sol%v_dn = this%vr_jm_fn(1) - this%Raf * ( this%qr_jm_fn(1) - flux(1:this%jms) )
    this%sol%v_up = this%vr_jm_fn(this%nd)
    
    this%sol%v_dn(1) = czero
    this%sol%v_up(1) = czero
    
    this%sol%u_dn = this%sol%u_dn + this%sol%v_dn * this%dt
    this%sol%u_up = this%sol%u_up + this%sol%v_up * this%dt
    
    do concurrent ( ijm = 2:this%jms )
      this%sol%t_dn(ijm) = this%sol%u_dn(ijm) - this%Vdelta_fn(1      ,ijm)
      this%sol%t_up(ijm) = this%sol%u_up(ijm) - this%Vdelta_fn(this%nd,ijm)
    end do
    
  end subroutine EE_topo_iceCrust_sub
  
end module IceCrustMod