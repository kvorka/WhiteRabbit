submodule (IceCrustMod) EE_iceCrust
  implicit none; contains
  
  module subroutine EE_iceCrust_sub(this, flux_bnd)
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
      call this%set_visc_sub()
      call this%set_lambda_sub()
      call this%set_cp_sub()
      call this%set_alpha_sub()
    
    call this%EE_mech_sub(flux)
    
    this%sol%v_dn(1) = czero
    this%sol%v_up(1) = czero
    
    do concurrent ( ijm = 2:this%jms )
      this%sol%v_dn(ijm) = this%vr_r_fn(1,ijm) - this%Raf * ( this%qr_r_fn(1,ijm) - flux(ijm) )
      this%sol%v_up(ijm) = this%vr_r_fn(this%nd,ijm)
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
  
  module subroutine EE_temp_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(inout) :: flux(:)
    integer                          :: ir, ijm
    
    !! At first, solve for degree zero in order to find the new heat flux
    ijm = 1
      ir = 1
        this%rtemp(ir,ijm) = cs4pi
      
      do concurrent ( ir = 2:this%nd )
        this%rtemp(ir,ijm) = this%htide_rr_fn(ir,ijm) + this%ntemp(ijm,ir)
      end do
      
      ir = this%nd+1
        this%rtemp(ir,ijm) = czero
      
      call this%solve_temp_sub( ijmstart=ijm, ijmend=ijm, ijmstep=1, rematrix=.true., matxsol=.true. )
    
    !! Update the heat flux
    flux = flux * c2r_fn( -this%q_r_fn(1,1,1) / s4pi )
    
    !! Solve for other degrees
    !$omp parallel do private (ir)
    do ijm = 2, this%jms
      ir = 1
        this%rtemp(1,ijm) = -( this%sol%u_dn(ijm) + ( this%vr_r_fn(1,ijm) + this%Raf * flux(ijm) ) * this%dt +        &
                             & this%Cl / ( c2r_fn( this%dT_dr_r_fn(ir,1) ) / s4pi - this%Cl ) * this%Vdelta_fn(1,ijm) )
      
      do concurrent ( ir = 2:this%nd )
        this%rtemp(ir,ijm) = this%htide_rr_fn(ir,ijm) + this%ntemp(ijm,ir)
      end do
      
      ir = this%nd+1
        this%rtemp(this%nd+1,ijm) = -( this%sol%u_up(ijm) + this%vr_r_fn(this%nd,ijm) * this%dt )
    end do
    !$omp end parallel do
    
    call this%solve_temp_sub( ijmstart=2, ijmend=this%jms, ijmstep=1, rematrix=.true., matxsol=.true. )
    
  end subroutine EE_temp_iceCrust_sub
  
  module subroutine EE_mech_iceCrust_sub(this, flux)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux(:)
    integer                          :: ir, ij, ijm
    complex(kind=dbl)                :: buoy
    
    !$omp parallel do private (ir,ij,buoy)
    do ijm = 2, this%jms
      ij = this%j_indx(ijm)

      ir = 1
        this%rsph1(ir,ijm) = -( this%sol%u_dn(ijm) - this%Vdelta_fn(1,ijm) -            &
                              & this%Raf * (this%qr_r_fn(ir,ijm) - flux(ijm)) * this%dt )
        this%rsph2(ir,ijm) = czero
      
      do ir = 2, this%nd
        buoy = this%Ra * this%alpha_rr_fn(ir) * this%gravity%g_fn( this%rad_grid%rr(ir) ) * this%temp_rr_fn(ir,ijm)
        
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
  
end submodule EE_iceCrust