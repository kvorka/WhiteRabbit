module IceCrustMod
  use IceMod
  use IceTidesMod
  implicit none
  
  type, extends(T_ice), public :: T_iceCrust
    type(T_iceTides), private :: tides
    
    contains
    
    procedure, public, pass :: init_sub        => init_iceCrust_sub
    procedure, public, pass :: time_scheme_sub => iter_iceCrust_sub
    procedure, public, pass :: solve_sub       => solve_iceCrust_sub
    
    procedure, private, pass :: EE_sub        => EE_iceCrust_sub
    procedure, private, pass :: EE_temp_sub   => EE_temp_iceCrust_sub
    procedure, private, pass :: EE_mech_sub   => EE_mech_iceCrust_sub
    procedure, private, pass :: adjust_dt_sub => adjust_dt_iceCrust_sub
    procedure, private, pass :: relevant_criterion_fn
    
    procedure, private, pass :: II_stress_fn => II_stress_iceCrust_fn
    procedure, private, pass :: avrg_temp_fn => avrg_temp_iceCrust_fn
    
  end type T_iceCrust
  
  private :: init_iceCrust_sub
  private :: iter_iceCrust_sub
  private :: solve_iceCrust_sub
  
  private :: vypis_iceCrust_sub
  private :: EE_iceCrust_sub
  private :: EE_temp_iceCrust_sub
  private :: EE_mech_iceCrust_sub
  private :: adjust_dt_iceCrust_sub
  private :: relevant_criterion_fn
  
  private :: II_stress_iceCrust_fn
  private :: avrg_temp_iceCrust_fn
  
  contains
  
  subroutine init_iceCrust_sub(this)
    class(T_iceCrust), intent(inout) :: this
   
    call this%init_ice_sub(jmax_in = jmax_ice, rheol_in = 'viscel', n_iter = n_iter_ice)
    call this%lat_grid%init_vcvv_sub()
    
    this%cf = 1._dbl

    call this%init_eq_temp_sub( rhs=.true. , nl=.true.  )
    call this%init_eq_mech_sub( rhs=.true. , nl=.false. )
    
    call this%tides%init_sub()
    
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
    
    do n = 1, this%n_iter
      call this%EE_sub( flux_bnd )
    end do

    call vypis_iceCrust_sub(this)
    
  end subroutine iter_iceCrust_sub
  
  subroutine solve_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux(:)
    integer                          :: iter, ir
    complex(kind=dbl), allocatable   :: Temp1(:), Temp2(:)
    
    allocate( Temp1(this%nd+1), Temp2(this%nd+1) )
    
    !! Seek for conductive solution with zero rhs at first
    this%dt = huge(0._dbl)
      do
        Temp1 = this%sol%temp_i_fn(1)
        
        ir = 1
          this%rtemp(ir,1) = cs4pi
        
        do concurrent ( ir = 2:this%nd )
          this%rtemp(ir,1) = czero
        end do
        
        ir = this%nd+1
          this%rtemp(ir,1) = czero
        
        call this%solve_temp_sub( ijmstart=1, ijmend=1, ijmstep=1, rematrix=.true., matxsol=.false. )
        
        if ( maxval(abs(this%sol%temp_i_fn(1) - Temp1)/abs(Temp1)) < 1e-8 ) exit
      end do    
    
    do iter = 1, 5
      !! Find tidal heating for given temperature rhs and stress field
      this%dt = huge(0._dbl)
        do
          call this%tides%compute_sub( this%II_stress_fn(), this%avrg_temp_fn() )
          this%htide = this%tides%htide
          
          Temp2 = this%sol%temp_i_fn(1)
          
          do
            Temp1 = this%sol%temp_i_fn(1)
            
            ir = 1
              this%rtemp(ir,1) = cs4pi
            
            do concurrent ( ir = 2:this%nd )
              this%rtemp(ir,1) = this%htide_fn(ir,1) + this%ntemp(1,ir)
            end do
            
            ir = this%nd+1
              this%rtemp(ir,1) = czero
            
            call this%solve_temp_sub( ijmstart=1, ijmend=1, ijmstep=1, rematrix=.true., matxsol=.false. )
            
            if ( maxval(abs(this%sol%temp_i_fn(1) - Temp1)/abs(Temp1)) < 1e-8 ) exit
          end do
          
          if ( maxval(abs(this%sol%temp_i_fn(1) - Temp2)/abs(Temp2)) < 1e-6 ) exit
        end do
      
      !! Solve for hydrostatic state for given tidal heating
      call this%set_dt_sub() ; call this%sol%nulify_sub()
        do
          call this%EE_sub(flux)
          
          if ( this%relevant_criterion_fn() < 5e-4 ) then 
            exit
          else
            call this%adjust_dt_sub()
          end if
        end do
        
      !! Vypisova kontrola
      write(*,*) this%sol%u_up(4)
    end do
    
    call this%set_dt_sub()
    
  end subroutine solve_iceCrust_sub
  
    function II_stress_iceCrust_fn(this) result(II_stress)
      class(T_iceCrust), intent(in) :: this
      integer                       :: ir
      real(kind=dbl), allocatable   :: II_stress(:)
      
      allocate( II_stress(this%nd) )
      
      do ir = 1, this%nd
        II_stress(ir) = tnorm_fn( this%jmax, this%sol%deviatoric_stress_jml2_fn(ir) )
      end do
      
      II_stress = II_stress * (this%viscU * this%kappaU / this%D_ud**2) / sqrt(4*pi)
      
    end function II_stress_iceCrust_fn
  
    function avrg_temp_iceCrust_fn(this) result (avrg_temp)
      class(T_iceCrust), intent(in) :: this
      integer                       :: ir
      real(kind=dbl), allocatable   :: avrg_temp(:)
      
      allocate( avrg_temp(this%nd) )
      
      do ir = 1, this%nd
        avrg_temp(ir) = this%Tu + (this%Td-this%Tu) * c2r_fn( this%rad_grid%c(ir,-1) * this%sol%temp_fn(ir  ,1) + &
                                                            & this%rad_grid%c(ir,+1) * this%sol%temp_fn(ir+1,1)   ) / sqrt(4*pi)
      end do
      
    end function avrg_temp_iceCrust_fn
  
  subroutine EE_iceCrust_sub(this, flux_bnd)
    class(T_iceCrust),           intent(inout) :: this
    complex(kind=dbl), optional, intent(in)    :: flux_bnd(:)
    integer                                    :: ir, ijm
    complex(kind=dbl), allocatable             :: flux(:)
    
    this%t = this%t + this%dt
    
    !$omp parallel do
    do ir = 2, this%nd
      this%ntemp(:,ir) = -this%vgradT_fn(ir)
    end do
    !$omp end parallel do
    
    allocate( flux(this%jms) )
      if ( present(flux_bnd) ) then
        flux = flux_bnd(1:this%jms)
      else
        flux = czero
      end if
      
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
    
    deallocate( flux )
    
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
  
  subroutine adjust_dt_iceCrust_sub(this)
    class(T_iceCrust), intent(inout) :: this
    
    if (this%dt < 0.2_dbl) then
      if ( this%relevant_criterion_fn() < 1.0d-2 ) then
        this%dt = 5 * this%dt
      end if
    else
      this%dt = 0.98_dbl
    end if
    
  end subroutine adjust_dt_iceCrust_sub
  
    real(kind=dbl) function relevant_criterion_fn(this)
      class(T_iceCrust), intent(inout) :: this
      real(kind=dbl)                   :: bottom, surface
      
      bottom = max( abs(this%sol%v_dn( 4) / this%sol%u_dn( 4)) , &
                  & abs(this%sol%v_dn( 6) / this%sol%u_dn( 6)) , &
                  & abs(this%sol%v_dn(11) / this%sol%u_dn(11))   ) * this%dt
      
      surface = max( abs(this%sol%v_up( 4)) / abs(this%sol%u_up( 4)) , &
                   & abs(this%sol%v_up( 6)) / abs(this%sol%u_up( 6)) , &
                   & abs(this%sol%v_up(11)) / abs(this%sol%u_up(11))   ) * this%dt
                   
      relevant_criterion_fn = surface
      
    end function relevant_criterion_fn
  
end module IceCrustMod