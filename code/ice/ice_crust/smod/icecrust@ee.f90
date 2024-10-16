submodule (icecrust) ee
  implicit none; contains
  
  module subroutine EE_iceCrust_sub(this, flux_bnd)
    class(T_iceCrust),           intent(inout) :: this
    complex(kind=dbl), optional, intent(in)    :: flux_bnd(:)
    integer                                    :: ijm
    complex(kind=dbl), allocatable             :: flux(:)
    
    this%t = this%t + this%dt
    
    call this%mvgradT_cpdivq_sub()
    
    allocate( flux(this%jms) )
      if ( present(flux_bnd) ) then
        do concurrent ( ijm = 1:this%jms )
          flux(ijm) = flux_bnd(ijm)
        end do
      else
        do concurrent ( ijm = 1:this%jms )
          flux(ijm) = czero
        end do
      end if
      
    call this%EE_temp_sub(flux)
      call this%set_visc_sub()
      call this%set_lambda_sub()
      call this%set_cp_sub()
      call this%set_alpha_sub()
    
    call this%EE_mech_sub(flux)
    
    this%bnd%v_dn(1) = czero
    this%bnd%v_up(1) = czero
    
    !$omp parallel do
    do ijm = 2, this%jms
      !Bottom bnd
      this%bnd%v_dn(ijm) = this%vr_r_fn(1,ijm) - this%Raf * ( this%qr_r_fn(1,ijm) - flux(ijm) )
      this%bnd%u_dn(ijm) = this%bnd%u_dn(ijm) + this%bnd%v_dn(ijm) * this%dt
      
      !Upper bnd
      this%bnd%v_up(ijm) = this%vr_r_fn(this%nd,ijm)
      this%bnd%u_up(ijm) = this%bnd%u_up(ijm) + this%bnd%v_up(ijm) * this%dt
    end do
    !$omp end parallel do
    
    !$omp parallel do
    do ijm = 2, this%jms
      this%bnd%t_dn(ijm) = this%bnd%u_dn(ijm) - this%Vdelta_fn(1      ,ijm)
      this%bnd%t_up(ijm) = this%bnd%u_up(ijm) - this%Vdelta_fn(this%nd,ijm)
    end do
    !$omp end parallel do
    
    deallocate( flux )
    
  end subroutine EE_iceCrust_sub
  
  module subroutine EE_temp_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(inout) :: flux(:)
    integer                          :: ir, ijm
    real(kind=dbl)                   :: gradTCoeff
    
    !! At first, solve for degree zero in order to find the new heat flux
    ijm = 1
      ir = 1
        this%rtemp(ir,ijm) = cs4pi
      
      do concurrent ( ir = 2:this%nd )
        this%rtemp(ir,ijm) = this%htide_rr_fn(ir,ijm) + this%ntemp(ijm,ir)
      end do
      
      ir = this%nd+1
        this%rtemp(ir,ijm) = this%bnd%temp_up(1)
      
      call this%solve_temp_sub( ijmstart=ijm, ijmend=ijm, ijmstep=1, rematrix=.true., matxsol=.true. )
    
    !! Update the heat flux
    flux = flux * c2r_fn( this%qr_r_fn(1,1) / s4pi )
    
    !! Solve for other degrees
    !$omp parallel do private (ir, gradTCoeff)
    do ijm = 2, this%jms
      ir = 1
        gradTCoeff = c2r_fn( this%dT_dr_r_fn(1,1) ) / s4pi - this%Cl
        this%rtemp(1,ijm) = -( this%bnd%u_dn(ijm) + ( this%vr_r_fn(1,ijm) + this%Raf * flux(ijm) ) * this%dt + &
                             & this%Cl / gradTCoeff * this%Vdelta_fn(1,ijm) )
      
      do concurrent ( ir = 2:this%nd )
        this%rtemp(ir,ijm) = this%htide_rr_fn(ir,ijm) + this%ntemp(ijm,ir)
      end do
      
      ir = this%nd+1
        gradTCoeff = c2r_fn( this%dT_dr_r_fn(this%nd,1) ) / s4pi
        this%rtemp(this%nd+1,ijm) = this%bnd%temp_up(ijm) / gradTCoeff - &
                                  & ( this%bnd%u_up(ijm) + this%vr_r_fn(this%nd,ijm) * this%dt )
    end do
    !$omp end parallel do
    
    call this%solve_temp_sub( ijmstart=2, ijmend=this%jms, ijmstep=1, rematrix=.true., matxsol=.true. )
    
  end subroutine EE_temp_iceCrust_sub
  
  module subroutine EE_mech_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux(:)
    integer                          :: ir, ij, ijm
    
    !$omp parallel do private (ir,ij)
    do ijm = 2, this%jms
      ij = this%j_indx(ijm)

      ir = 1
        this%rsph1(ir,ijm) = -( this%bnd%u_dn(ijm) - this%Vdelta_fn(1,ijm) -            &
                              & this%Raf * (this%qr_r_fn(ir,ijm) - flux(ijm)) * this%dt )
        this%rsph2(ir,ijm) = czero
      
      do ir = 2, this%nd
        this%rsph1(ir,ijm) = this%buoy_rr_fn(ir,-1,ijm,-1)
        this%rsph2(ir,ijm) = this%buoy_rr_fn(ir,+1,ijm,-1)
      end do
      
      ir = this%nd+1
        this%rsph1(ir,ijm) = czero
        this%rsph2(ir,ijm) = -( this%bnd%u_up(ijm) - this%Vdelta_fn(this%nd,ijm) )
    end do
    !$omp end parallel do

    call this%solve_mech_sub( ijmstart=2, ijmend=this%jms, ijmstep=1, rematrix=.true., matxsol=.true. )
    
  end subroutine EE_mech_iceCrust_sub
  
end submodule ee