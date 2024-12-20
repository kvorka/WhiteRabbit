submodule (physicalobject) equations_mech
  implicit none ; contains
  
  module procedure init_eq_mech_sub
    
    call this%sol%init_smech_sub()
    call this%mat%init_mmech_sub()

    allocate( this%rsph1(this%nd+1,this%jms) ) ; this%rsph1 = czero
    allocate( this%rsph2(this%nd+1,this%jms) ) ; this%rsph2 = czero
    
  end procedure init_eq_mech_sub
  
  module procedure prepare_mat_mech_sub
    integer :: ij
    
    !$omp parallel do
    do ij = ijstart, ijend
      call this%mat%mech(ij)%fill_sub( this%mat_mech_fn(j_in=ij, a_in=this%cf  ), &
                                     & this%mat_mech_fn(j_in=ij, a_in=this%cf-1)  )
    end do
    !$omp end parallel do
    
  end procedure prepare_mat_mech_sub
  
  module procedure solve_mech_sub
    integer :: ij, ijm, ir, is
    
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
            this%sol%mech(is+3,ijm) = this%Ramu * this%sol%mech(is+2,ijm) / this%dt
            this%sol%mech(is+2,ijm) = czero
            this%sol%mech(is+4,ijm) = this%Ramu * this%sol%mech(is+4,ijm) / this%dt
            this%sol%mech(is+5,ijm) = this%Ramu * this%sol%mech(is+5,ijm) / this%dt
          end do
            
          ir = this%nd+1
            this%sol%mech(6*this%nd+1,ijm) = this%rsph1(ir,ijm)
            this%sol%mech(6*this%nd+2,ijm) = this%rsph2(ir,ijm)
            
          call this%mat%mech(ij)%luSolve_sub( this%sol%mech(:,ijm) )
        end do
        !$omp end parallel do
      end select
      
  end procedure solve_mech_sub
  
end submodule equations_mech