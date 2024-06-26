module IceCrustMod
  use IceMod
  use IceTidesMod
  implicit none
  
  type, extends(T_ice), public :: T_iceCrust
    type(T_iceTides), private :: tides
    
    contains
    
    procedure :: init_sub        => init_iceCrust_sub
    procedure :: time_scheme_sub => iter_iceCrust_sub
    procedure :: solve_sub       => solve_iceCrust_sub
    
    procedure, private :: EE_sub        => EE_iceCrust_sub
    procedure, private :: EE_temp_sub   => EE_temp_iceCrust_sub
    procedure, private :: EE_mech_sub   => EE_mech_iceCrust_sub
    procedure, private :: II_stress_fn => II_stress_iceCrust_fn
    procedure, private :: avrg_temp_fn => avrg_temp_iceCrust_fn
    
  end type T_iceCrust
  
  private :: init_iceCrust_sub, iter_iceCrust_sub, solve_iceCrust_sub, vypis_iceCrust_sub, EE_iceCrust_sub, EE_temp_iceCrust_sub, &
           & EE_mech_iceCrust_sub, II_stress_iceCrust_fn, avrg_temp_iceCrust_fn
  
  contains
  
  subroutine init_iceCrust_sub(this)
    class(T_iceCrust), intent(inout) :: this
   
    call this%init_ice_sub(jmax_in = jmax_ice, rheol_in = 'viscel', n_iter = n_iter_ice)
    
    this%cf = one

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
    class(T_iceCrust),           intent(inout) :: this
    complex(kind=dbl), optional, intent(in)    :: flux(:)
    integer                                    :: iter, ir
    complex(kind=dbl), allocatable             :: Temp1(:), Temp2(:), u_up1(:)
    
    allocate( Temp1(this%nd+1), Temp2(this%nd+1), u_up1(this%jms) )
    
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
    
    do
      !! Save latest value of shape
      u_up1(:) = this%sol%u_up(:)

      !! Find tidal heating for given temperature rhs and stress field
      this%dt = huge(zero)
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
      
      !! Solve for given tidal heating
      call this%set_dt_sub() ; call this%sol%nulify_sub()
        do
          if ( present(flux) ) then
            call this%EE_sub(flux_bnd=flux)
          else
            call this%EE_sub()
          end if
          
          if ( max( abs( this%sol%v_up( 4) * this%dt / this%sol%u_up( 4) ) , &
                  & abs( this%sol%v_up( 6) * this%dt / this%sol%u_up( 6) ) , &
                  & abs( this%sol%v_up(11) * this%dt / this%sol%u_up(11) ) , &   
                  & abs( this%sol%v_up(13) * this%dt / this%sol%u_up(13) ) , & 
                  & abs( this%sol%v_up(15) * this%dt / this%sol%u_up(15) )   ) < 1e-4 ) then 
            exit
          else if ( this%dt < 0.2_dbl ) then
              if ( max( abs( this%sol%v_up( 4) * this%dt / this%sol%u_up( 4) ) , &
                      & abs( this%sol%v_up( 6) * this%dt / this%sol%u_up( 6) ) , &
                      & abs( this%sol%v_up(11) * this%dt / this%sol%u_up(11) ) , &   
                      & abs( this%sol%v_up(13) * this%dt / this%sol%u_up(13) ) , & 
                      & abs( this%sol%v_up(15) * this%dt / this%sol%u_up(15) )   ) < 1e-3 ) this%dt = 2 * this%dt
          else
            this%dt = 0.48_dbl
          end if
        end do
      
      !! Stopping criterion
      if ( max ( abs( (u_up1( 4)-this%sol%u_up( 4))/this%sol%u_up( 4) )  ,          &
         &       abs( (u_up1( 6)-this%sol%u_up( 6))/this%sol%u_up( 6) )  ,          &
         &       abs( (u_up1(11)-this%sol%u_up(11))/this%sol%u_up(11) )  ,          &
         &       abs( (u_up1(13)-this%sol%u_up(13))/this%sol%u_up(13) )  ,          &
         &       abs( (u_up1(15)-this%sol%u_up(15))/this%sol%u_up(15) )  ) < 1e-4 ) exit
    end do
    
    deallocate( Temp1, Temp2, u_up1 ) ; call this%set_dt_sub()
    
  end subroutine solve_iceCrust_sub
  
    function II_stress_iceCrust_fn(this) result(II_stress)
      class(T_iceCrust), intent(in)  :: this
      integer                        :: ir
      real(kind=dbl)                 :: facstress
      real(kind=dbl),    allocatable :: II_stress(:)
      complex(kind=dbl), allocatable :: devtens(:)
      
      facstress = (this%viscU * this%kappaU / this%D_ud**2) / s4pi
      
      allocate( II_stress(this%nd), devtens(this%jmt) )
      
      do ir = 1, this%nd
        devtens       = this%sol%deviatoric_stress_jml2_fn(ir)
        II_stress(ir) = sqrt( tensproduct_fn(this%jmax, devtens, devtens) ) * facstress
      end do
      
      deallocate( devtens )
      
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
      call this%mvgradT_sub(ir, this%ntemp(:,ir))
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
    
    this%sol%v_dn(1) = czero
    this%sol%v_up(1) = czero
    
    do concurrent ( ijm = 2:this%jms )
      this%sol%v_dn(ijm) = this%vr_fn(1,ijm) - this%Raf * ( this%qr_fn(1,ijm) - flux(ijm) )
      this%sol%v_up(ijm) = this%vr_fn(this%nd,ijm)
    end do
    
    do concurrent ( ijm = 1:this%jms )
      this%sol%u_dn(ijm) = this%sol%u_dn(ijm) + this%sol%v_dn(ijm) * this%dt
      this%sol%u_up(ijm) = this%sol%u_up(ijm) + this%sol%v_up(ijm) * this%dt
    end do
    
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
    real(kind=dbl)                   :: rgrad_T
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
    
    !$omp parallel do private (ir, rgrad_T)
    do ijm = 2, this%jms
      rgrad_T = - ( c2r_fn( -this%sol%flux_fn(1,1,1) / sqrt(4*pi) ) / this%lambda_fn(1) )
      
      ir = 1
        this%rtemp(1,ijm) = -( this%sol%u_dn(ijm) + ( this%vr_fn(1,ijm) + this%Raf * flux(ijm) ) * this%dt + &
                             & this%Cl / ( rgrad_T - this%Cl ) * this%Vdelta_fn(1,ijm)                       )
      
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
    complex(kind=dbl)                :: buoy
    
    !$omp parallel do private (ir,ij,buoy)
    do ijm = 2, this%jms
      ij = this%j_indx(ijm)

      ir = 1
        this%rsph1(ir,ijm) = -( this%sol%u_dn(ijm) - this%Vdelta_fn(1,ijm) - this%Raf * (this%qr_fn(ir,ijm) - flux(ijm)) * this%dt )
        this%rsph2(ir,ijm) = czero
      
      do ir = 2, this%nd
        buoy = this%Ra * this%alpha_fn(ir) * this%gravity%g_fn( this%rad_grid%rr(ir) ) * this%sol%temp_fn(ir,ijm)
        
        this%rsph1(ir,ijm) = -sqrt((ij  )/(2*ij+one)) * buoy
        this%rsph2(ir,ijm) = +sqrt((ij+1)/(2*ij+one)) * buoy
      end do
      
      ir = this%nd+1
        this%rsph1(ir,ijm) = czero
        this%rsph2(ir,ijm) = -( this%sol%u_up(ijm) - this%Vdelta_fn(this%nd,ijm) )
    end do
    !$omp end parallel do

    call this%solve_mech_sub( ijmstart=2, ijmend=this%jms, ijmstep=1, rematrix=.true., matxsol=.true. )
    
  end subroutine EE_mech_iceCrust_sub
  
end module IceCrustMod