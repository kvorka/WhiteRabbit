submodule (PhysicalObject) Equations_mech
  implicit none
  contains

  subroutine init_eq_mech_sub(this, rhs, nl)
    class(T_physicalObject), intent(inout) :: this
    logical,                 intent(in)    :: rhs, nl
    
    call this%sol%init_smech_sub()
    call this%mat%init_mmech_sub()

    if (rhs) then
      allocate( this%rsph1(this%nd+1,this%jms) ) ; this%rsph1 = czero
      allocate( this%rsph2(this%nd+1,this%jms) ) ; this%rsph2 = czero
    end if

    if (nl) then
      allocate( this%nsph1(this%jms,2:this%nd) ) ; this%nsph1 = czero
      allocate( this%nsph2(this%jms,2:this%nd) ) ; this%nsph2 = czero
    end if
    
  end subroutine init_eq_mech_sub
  
  subroutine prepare_mat_mech_sub(this, ijstart, ijend)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: ijstart, ijend
    integer                                :: ij
    
    !$omp parallel do
    do ij = ijstart, ijend
      call this%mat%mech(ij)%fill_sub( this%matica_mech_fn(j_in=ij, a_in=this%cf  ), &
                                     & this%matica_mech_fn(j_in=ij, a_in=this%cf-1)  )
    end do
    !$omp end parallel do
    
  end subroutine prepare_mat_mech_sub
  
  subroutine solve_mech_sub(this, ijmstart, ijmend, ijmstep, rematrix, matxsol)
    class(T_physicalObject), intent(inout)        :: this
    integer,                 intent(in)           :: ijmstart, ijmend, ijmstep
    logical,                 intent(in)           :: rematrix, matxsol
    integer                                       :: ij, ijm, ir, is
    real(kind=dbl), allocatable                   :: visc(:)
    
    if (rematrix) call this%prepare_mat_mech_sub( this%j_indx(ijmstart), this%j_indx(ijmend) )
    
    select case (this%rheology)
      case ('viscos')
        !$omp parallel do private (ir,is,ij)
        do ijm = ijmstart, ijmend, ijmstep
          ij = this%j_indx(ijm)
          
          if ( matxsol ) then
            do concurrent ( ir=2:this%nd )
              this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + this%mat%mech(ij)%multipl_fn(6*(ir-1)+1,this%sol%mech(:,ijm))
              this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + this%mat%mech(ij)%multipl_fn(6*(ir-1)+2,this%sol%mech(:,ijm))
            end do
          end if
          
          do concurrent ( ir=1:this%nd )
            is = 6*(ir-1)+1
            
            this%sol%mech(is  ,ijm) = this%rsph1(ir,ijm)
            this%sol%mech(is+1,ijm) = this%rsph2(ir,ijm)
            this%sol%mech(is+2,ijm) = czero
            this%sol%mech(is+3,ijm) = czero
            this%sol%mech(is+4,ijm) = czero
            this%sol%mech(is+5,ijm) = czero
          end do
            
          ir = this%nd+1
            this%sol%mech(6*this%nd+1,ijm) = this%rsph1(ir,ijm)
            this%sol%mech(6*this%nd+2,ijm) = this%rsph2(ir,ijm)
            
          call this%mat%mech(ij)%luSolve_sub( this%sol%mech(:,ijm) )
        end do
        !$omp end parallel do
      
      case ('viscel')
        allocate( visc(this%nd) )
          do concurrent ( ir=1:this%nd )
            visc(ir) = this%visc_fn(ir)
          end do
        
        !$omp parallel do private (ir,is,ij)
        do ijm = ijmstart, ijmend, ijmstep
          ij = this%j_indx(ijm)
          
          if ( matxsol ) then
            do concurrent ( ir=2:this%nd )
              this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + this%mat%mech(ij)%multipl_fn(6*(ir-1)+1,this%sol%mech(:,ijm))
              this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + this%mat%mech(ij)%multipl_fn(6*(ir-1)+2,this%sol%mech(:,ijm))
            end do
          end if
          
          do concurrent ( ir=1:this%nd )
            is = 6*(ir-1)+1
            
            this%sol%mech(is  ,ijm) = this%rsph1(ir,ijm)
            this%sol%mech(is+1,ijm) = this%rsph2(ir,ijm)
            this%sol%mech(is+3,ijm) = this%Ramu * visc(ir) * this%sol%mech(is+2,ijm) / this%dt
            this%sol%mech(is+2,ijm) = czero
            this%sol%mech(is+4,ijm) = this%Ramu * visc(ir) * this%sol%mech(is+4,ijm) / this%dt
            this%sol%mech(is+5,ijm) = this%Ramu * visc(ir) * this%sol%mech(is+5,ijm) / this%dt
          end do
            
          ir = this%nd+1
            this%sol%mech(6*this%nd+1,ijm) = this%rsph1(ir,ijm)
            this%sol%mech(6*this%nd+2,ijm) = this%rsph2(ir,ijm)
            
          call this%mat%mech(ij)%luSolve_sub( this%sol%mech(:,ijm) )
        end do
        !$omp end parallel do
        
        deallocate( visc )
      end select
      
  end subroutine solve_mech_sub

end submodule Equations_mech