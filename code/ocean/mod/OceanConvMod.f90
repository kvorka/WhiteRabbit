module OceanConvMod
  use OceanMod
  implicit none
  
  type, extends(T_ocean), public :: T_oceanConv
    contains
    
    procedure, public, pass :: init_sub        => init_oceanConv_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanConv_sub
  end type T_oceanConv
  
  contains
  
  subroutine init_oceanConv_sub(this)
    class(T_oceanConv), intent(inout) :: this
    
    call this%init_ocean_sub()
    
    call this%init_eq_temp_sub( rhs=.true. , nl=.true. )
    call this%init_eq_torr_sub( rhs=.true. , nl=.true. )
    call this%init_eq_mech_sub( rhs=.true. , nl=.true. )
    
    call this%init_state_sub()
    
  end subroutine init_oceanConv_sub
  
  subroutine time_scheme_oceanConv_sub(this)
    class(T_oceanConv), intent(inout) :: this
    integer                           :: ir, ijm
    real(kind=dbl)                    :: dt
    logical                           :: changed_dt
    
    changed_dt = .false. ; dt = this%velc_crit_fn()
    
    if ( dt < this%dt) then
      dt = 0.50_dbl * dt
      
      this%ab = 1 + dt / ( 2 * this%dt )
      this%dt = dt
      
      changed_dt = .true.
    else
      this%ab    = 1.5_dbl
      changed_dt = .false.
    end if
    
    ijm = 1 ; ir = 1
      this%rtemp(ir,ijm) = cs4pi
    
    !$omp parallel
    !$omp do collapse (2)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtemp(ir,ijm) = (1-this%ab) * this%ntemp(ijm,ir)
        this%rtorr(ir,ijm) = (1-this%ab) * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = (1-this%ab) * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = (1-this%ab) * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end do
    
    !$omp do
    do ir = 2, this%nd
      call this%fullnl_sub(ir)
    end do
    !$omp end do
    
    !$omp do collapse (2)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtemp(ir,ijm) = this%rtemp(ir,ijm) + this%ab * this%ntemp(ijm,ir)
        this%rtorr(ir,ijm) = this%rtorr(ir,ijm) + this%ab * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + this%ab * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + this%ab * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end do
    !$omp end parallel
    
    call this%solve_temp_sub( ijmstart=1 , ijmend=this%jms, ijmstep=1, rematrix=changed_dt, matxsol=.true. )
    call this%solve_torr_sub( ijmstart=2 , ijmend=this%jms, ijmstep=1, rematrix=changed_dt, matxsol=.true. )
    call this%solve_mech_sub( ijmstart=2 , ijmend=this%jms, ijmstep=1, rematrix=changed_dt, matxsol=.true. )
    
    if (this%mechanic_bnd == 'frees') call this%global_rotation_sub()
    
  end subroutine time_scheme_oceanConv_sub
  
end module OceanConvMod
