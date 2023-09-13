module IceCrustMod
  use IceMod
  use IceTidesMod
  implicit none
  
  type, extends(T_ice), public :: T_iceCrust
    contains
    
    procedure, public, pass :: init_sub        => init_iceCrust_sub
    procedure, public, pass :: time_scheme_sub => iter_iceCrust_sub
    
    procedure, private, pass :: EE_sub        => EE_iceCrust_sub
    procedure, private, pass :: EE_temp_sub   => EE_temp_iceCrust_sub
    procedure, private, pass :: EE_mech_sub   => EE_mech_iceCrust_sub
    
  end type T_iceCrust
  
  private :: init_iceCrust_sub
  private :: vypis_iceCrust_sub
  
  private :: EE_iceCrust_sub
  private :: EE_temp_iceCrust_sub
  private :: EE_mech_iceCrust_sub
  
  contains
  
  subroutine init_iceCrust_sub(this, notides)
    class(T_iceCrust), intent(inout) :: this
    logical, optional, intent(in)    :: notides
    type(T_iceTides)                 :: tides
    integer                          :: i, iter
    real(kind=dbl),    allocatable   :: stress_i(:)
    complex(kind=dbl), allocatable   :: flux_bnd(:)
   
    call this%init_ice_sub(jmax_in = jmax_ice, rheol_in = 'viscel', n_iter = n_iter_ice)
    call this%lat_grid%init_vcvv_sub()
    
    this%cf = 1._dbl

    call this%init_eq_temp_sub( rhs=.true. , nl=.true.  )
    call this%init_eq_mech_sub( rhs=.true. , nl=.false. )
    
    call tides%init_sub()
    allocate( flux_bnd(this%jms) ); flux_bnd = czero
    allocate( stress_i(this%nd)  ); stress_i   = 0._dbl
    
      call tides%compute_sub()
      
      this%sol%temp(1:3*this%nd+1:3,1) = tides%sol%temp(1:3*this%nd+1:3,1)
      this%htide                       = tides%htide
      call this%set_dt_sub()
      
      iter = 0
        do
          call this%EE_sub(flux_bnd)
          write(*,*) this%dt, c2r_fn(this%sol%u_up(4)) * this%D_ud
          
          if ( abs(this%sol%v_up(4) / this%sol%u_up(4)) < 1e-10 ) then 
            exit
          else if ( abs(this%sol%v_up(4) * this%dt / this%sol%u_up(4)) < 1e-2 ) then
            this%dt = 5 * this%dt
          end if
        end do
      
      do iter = 1, 5
        do i = 1, this%nd
          stress_i(i) = tnorm_fn( this%jmax, this%sol%deviatoric_stress_jml2_fn(i) ) / sqrt(4*pi)
        end do
        
        call tides%compute_sub( stress_prof_i = (this%viscU * this%kappaU / this%D_ud**2) * stress_i )
        call this%sol%nulify_sub()
        
        this%sol%temp(1:3*this%nd+1:3,1) = tides%sol%temp(1:3*this%nd+1:3,1)
        this%htide                       = tides%htide
        call this%set_dt_sub()
        
        do
          call this%EE_sub(flux_bnd)
          write(*,*) this%dt, c2r_fn(this%sol%u_up(4)) * this%D_ud
          
          if ( abs(this%sol%v_up(4) / this%sol%u_up(4)) < 1e-10 ) then 
            exit
          else if ( abs(this%sol%v_up(4) * this%dt / this%sol%u_up(4)) < 1e-3 ) then
            this%dt = 5 * this%dt
          end if
        end do
      end do
      
    call tides%deallocate_sub()
    deallocate( flux_bnd, stress_i )
    
    write(*,*)
    write(*,*) c2r_fn(this%sol%u_up(4)) * this%D_ud , this%sol%u_up(6) * this%D_ud , c2r_fn(this%sol%u_up(11)) * this%D_ud
    
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
    integer                          :: ir, ijm
    
    this%t = this%t + this%dt
    
    !$omp parallel do
    do ir = 2, this%nd
      this%ntemp(:,ir) = -this%vgradT_fn(ir)
    end do
    !$omp end parallel do
    
    call this%EE_temp_sub(flux)
    call this%EE_mech_sub(flux)
    
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
    
  end subroutine EE_iceCrust_sub
  
  subroutine EE_temp_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(inout) :: flux(:)
    integer                          :: ir, ijm
    complex(kind=dbl), allocatable   :: Temp(:), Temp1(:)
    
    ijm = 1
      allocate( Temp(this%nd+1), Temp1(this%nd+1) ); Temp = this%sol%temp_i_fn(ijm)
        do
          Temp1 = this%sol%temp_i_fn(ijm)
          
          ir = 1
            this%rtemp(ir,ijm) = cs4pi
          
          do concurrent ( ir = 2:this%nd )
            this%rtemp(ir,ijm) = Temp(ir) / this%dt + this%htide_fn(ir,ijm) + this%ntemp(ijm,ir)
          end do
          
          ir = this%nd+1
            this%rtemp(ir,ijm) = czero
          
          call this%solve_temp_sub( ijmstart=ijm, ijmend=ijm, ijmstep=1, rematrix=.true., matxsol=.false. )
          
          if ( maxval(abs(this%sol%temp_i_fn(ijm) - Temp1)/abs(Temp1)) < 1e-5 ) exit       
        end do
      deallocate( Temp, Temp1 )
    
      flux = flux * c2r_fn( -this%sol%flux_fn(1,1,1) / sqrt(4*pi) )
    
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
    
    call this%solve_temp_sub( ijmstart=2, ijmend=this%jms, ijmstep=1, rematrix=.true., matxsol=.true. )
    
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

    call this%solve_mech_sub( ijmstart=2, ijmend=this%jms, ijmstep=1, rematrix=.true., matxsol=.true. )
    
  end subroutine EE_mech_iceCrust_sub
  
end module IceCrustMod