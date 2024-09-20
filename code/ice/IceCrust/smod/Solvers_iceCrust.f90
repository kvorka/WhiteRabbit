submodule (IceCrustMod) Solvers_iceCrust
  implicit none; contains
  
  module subroutine iter_iceCrust_sub(this, flux_bnd)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux_bnd(:)
    integer                          :: n
    
    do n = 1, this%n_iter
      call this%EE_sub( flux_bnd )
    end do

    call vypis_iceCrust_sub(this)
    
  end subroutine iter_iceCrust_sub
  
  module subroutine solve_iceCrust_sub(this, flux)
    class(T_iceCrust),           intent(inout) :: this
    complex(kind=dbl), optional, intent(in)    :: flux(:)
    integer                                    :: ir
    complex(kind=dbl), allocatable             :: u_up1(:)
    
    !! Seek for conductive solution with zero rhs at first
    call this%solve_conduction_sub()
    
    allocate( u_up1(this%jms) )
    
    do
      !! Save latest value of shape
      u_up1(:) = this%bnd%u_up(:)
      
      !! Find tidal heating for given viscosity field
      call this%tides%compute_sub( this%mparams%visc )
      call this%tides%htide_ir_ijm_sub( this%htide )
      
      !! Solve for given tidal heating
      call this%sol%nulify_sub()
      
      do
        if ( present(flux) ) then
          call this%EE_sub(flux_bnd=flux)
        else
          call this%EE_sub()
        end if
        
        if ( abs( this%bnd%v_up(4) * this%dt / this%bnd%u_up(4) ) < 1e-5 ) then 
          exit
        else if ( this%dt < 0.2_dbl ) then
            if ( abs( this%bnd%v_up(4) * this%dt / this%bnd%u_up(4) ) < 1e-4 ) this%dt = 2 * this%dt
        else
          this%dt = 0.48_dbl
        end if
      end do
      
      !! Stopping criterion
      write(*,*) abs( (u_up1(4)-this%bnd%u_up(4))/this%bnd%u_up(4) )
      if ( abs( (u_up1(4)-this%bnd%u_up(4))/this%bnd%u_up(4) ) < 1e-5 ) exit
    end do
    
    deallocate( u_up1 )
    
  end subroutine solve_iceCrust_sub
  
  module subroutine solve_conduction_iceCrust_sub(this)
    class(T_iceCrust), intent(inout) :: this
    integer                          :: ir, ijm
    complex(kind=dbl), allocatable   :: Temp1(:), Temp2(:)
    
    !!Seeking for conduction solution with zero rhs 
    !!does not require time stepping
    this%dt = huge(zero)
    
    !!Control array for mean temperature
    allocate( Temp1(this%nd+1), Temp2(this%nd+1) )
      
      do
        !! Save the mean temperature before the iteration
        call this%temp_irr_jm_sub(1, Temp1)
        
        !! Prepare the zero RHS equations
        !$omp parallel do private (ir)
        do ijm = 1, this%jms
          ir = 1
            if ( ijm > 1) then
              this%rtemp(1,ijm) = czero
            else
              this%rtemp(1,ijm) = cs4pi
            end if
            
          do concurrent ( ir = 2:this%nd+1 )
            this%rtemp(ir,ijm) = czero
          end do
        end do
        !$omp end parallel do
        
        !! Solve the conduction problem
        call this%solve_temp_sub( ijmstart=1, ijmend=this%jms, ijmstep=1, rematrix=.true., matxsol=.true. )
        
        !! Check the difference in mean temperature before and after the iteration
        call this%temp_irr_jm_sub(1, Temp2); if ( maxval(abs(Temp2-Temp1)/abs(Temp1)) < 1e-8 ) exit
        
        !! Update the parameters for the next iteration
        call this%set_lambda_sub()
        call this%set_cp_sub()
      end do
      
    deallocate( Temp1, Temp2 )
    
    !! Set back the time step
    call this%set_dt_sub()
    
    !! Set the material parameters
    call this%set_lambda_sub()
    call this%set_alpha_sub()
    call this%set_cp_sub()
    call this%set_visc_sub()
    
  end subroutine solve_conduction_iceCrust_sub
  
end submodule Solvers_iceCrust