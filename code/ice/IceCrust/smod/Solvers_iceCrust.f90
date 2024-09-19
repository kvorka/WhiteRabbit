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
    integer                                    :: iter, ir
    complex(kind=dbl), allocatable             :: Temp1(:), Temp2(:), Temp3(:), u_up1(:)
    
    !! Seek for conductive solution with zero rhs at first
    this%dt = huge(zero)
    allocate( Temp1(this%nd+1), Temp2(this%nd+1) )
      
      do
        call this%temp_irr_jm_sub(1, Temp1)
        
          !! Solve for given lambda
          call this%lambda_iceCrust_jm_sub()
          
          ir = 1
            this%rtemp(ir,1) = cs4pi
          
          do concurrent ( ir = 2:this%nd+1 )
            this%rtemp(ir,1) = czero
          end do
          
          call this%solve_temp_sub( ijmstart=1, ijmend=1, ijmstep=1, rematrix=.true., matxsol=.false. )
        
        call this%temp_irr_jm_sub(1, Temp2); if ( maxval(abs(Temp2 - Temp1)/abs(Temp1)) < 1e-8 ) exit
      end do
      
    deallocate( Temp1, Temp2 )
    
    allocate( Temp1(this%nd+1), Temp2(this%nd+1), Temp3(this%nd+1), u_up1(this%jms) )
    
    !! Start iterative solver from conductive solution
    do
      !! Save latest value of shape
      u_up1(:) = this%sol%u_up(:)

      !! Find tidal heating for given temperature rhs and stress field
      this%dt = huge(zero)
        do
          call this%temp_irr_jm_sub(1, Temp1)
          call this%tides%compute_sub( this%mparams%visc )
          this%htide = this%tides%htide
          
          do
            call this%temp_irr_jm_sub(1, Temp2)
            
            ir = 1
              this%rtemp(ir,1) = cs4pi
            
            do concurrent ( ir = 2:this%nd )
              this%rtemp(ir,1) = this%htide_fn(ir,1) + this%ntemp(1,ir)
            end do
            
            ir = this%nd+1
              this%rtemp(ir,1) = czero
            
            call this%solve_temp_sub( ijmstart=1, ijmend=1, ijmstep=1, rematrix=.true., matxsol=.false. )
            
            call this%temp_irr_jm_sub(1, Temp3) ; if ( maxval(abs(Temp3 - Temp2)/abs(Temp2)) < 1e-8 ) exit
          end do
          
          if ( maxval(abs(Temp2 - Temp1)/abs(Temp1)) < 1e-6 ) exit
        end do
      
      !! Solve for given tidal heating
      call this%set_dt_sub()
      call this%sol%nulify_sub()
      
        do
          if ( present(flux) ) then
            call this%EE_sub(flux_bnd=flux)
          else
            call this%EE_sub()
          end if
          
          if ( abs( this%sol%v_up(4) * this%dt / this%sol%u_up(4) ) < 1e-4 ) then 
            exit
          else if ( this%dt < 0.2_dbl ) then
              if ( abs( this%sol%v_up(4) * this%dt / this%sol%u_up(4) ) < 1e-3 ) this%dt = 2 * this%dt
          else
            this%dt = 0.48_dbl
          end if
        end do
      
      !! Stopping criterion
      write(*,*) abs( (u_up1(4)-this%sol%u_up(4))/this%sol%u_up(4) )
      if ( abs( (u_up1(4)-this%sol%u_up(4))/this%sol%u_up(4) ) < 1e-4 ) exit
    end do
    
    deallocate( Temp1, Temp2, u_up1 ) ; call this%set_dt_sub()
    
  end subroutine solve_iceCrust_sub
  
end submodule Solvers_iceCrust