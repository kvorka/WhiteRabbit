submodule (oceanconv) timescheme
  implicit none; contains
  
  module procedure time_scheme_oceanConv_sub
    integer :: ir, ijm
    
    !ir = 1 ; ijm = 1
      this%rtemp(1,1) = cs4pi
    
    !boundary conditions
    if ( this%init_thermal_bnd ) then
      
      select case ( this%thermal_bnd )
        case ('fluxd')
          do ijm = 2, this%jms
            this%rtemp(1,ijm) = this%bnd%flux_dn(ijm)
          end do
      
      end select
    end if
    
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
    !$omp end parallel
    
    call this%fullnl_sub()
    
    !$omp parallel do collapse (2)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtemp(ir,ijm) = this%rtemp(ir,ijm) + this%ab * this%ntemp(ijm,ir)
        this%rtorr(ir,ijm) = this%rtorr(ir,ijm) + this%ab * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + this%ab * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + this%ab * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end parallel do
    
    call this%solve_temp_sub( ijmstart=1 , ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.true. )
    call this%solve_torr_sub( ijmstart=2 , ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.true. )
    call this%solve_mech_sub( ijmstart=2 , ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.true. )
    
    if (this%mechanic_bnd == 'frees') call this%global_rotation_sub()
    
  end procedure time_scheme_oceanConv_sub
  
end submodule timescheme
