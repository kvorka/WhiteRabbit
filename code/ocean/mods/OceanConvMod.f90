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
      call this%lat_grid%init_vcsv_vcvv_vcvgv_sub()
    
    call this%init_eq_temp_sub(rhs=.true., nl=.true.)
    call this%init_eq_torr_sub(rhs=.true., nl=.true.)
    call this%init_eq_mech_sub(rhs=.true., nl=.true.)
    
    call this%init_state_sub()
    
  end subroutine init_oceanConv_sub
  
  subroutine time_scheme_oceanConv_sub(this, cf)
    class(T_oceanConv), intent(inout) :: this
    real(kind=dbl),     intent(in)    :: cf
    integer                           :: ir, ijm

    !$omp parallel
    !$omp do collapse (2)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtemp(ir,ijm) = (1-cf) * this%ntemp(ijm,ir)
        this%rtorr(ir,ijm) = (1-cf) * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = (1-cf) * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = (1-cf) * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end do

    !$omp do
    do ir = 2, this%nd
      call fullnl_sub( this, ir, this%ntemp(:,ir), &
                                 this%nsph1(:,ir), &
                                 this%ntorr(:,ir), &
                                 this%nsph2(:,ir)  )
    end do
    !$omp end do

    !$omp do collapse (2)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtemp(ir,ijm) = this%rtemp(ir,ijm) + cf * this%ntemp(ijm,ir)
        this%rtorr(ir,ijm) = this%rtorr(ir,ijm) + cf * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + cf * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + cf * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end do
    !$omp end parallel

    call this%solve_temp_sub()
    call this%solve_torr_sub()
    call this%solve_mech_sub()
    
    if (this%mechanic_bnd == 'frees') then
      call this%global_rotation_sub()
    end if
    
  end subroutine time_scheme_oceanConv_sub
  
end module OceanConvMod