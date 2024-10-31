submodule (physicalobject) equations_temp
  implicit none ; contains
  
  module procedure init_eq_temp_sub
    
    call this%sol%init_stemp_sub()
    call this%mat%init_mtemp_sub()
    
    allocate( this%rtemp(this%nd+1,this%jms) )
      this%rtemp = czero
    
  end procedure init_eq_temp_sub
  
  module procedure prepare_mat_temp_sub
    integer :: ij
    
    !$omp parallel do
    do ij = ijstart, ijend
      call this%mat%temp(ij)%fill_sub( this%mat_temp_fn(j_in=ij, a_in=this%cf),  this%mat_temp_fn(j_in=ij, a_in=this%cf-1) )
    end do
    !$omp end parallel do
    
  end procedure prepare_mat_temp_sub
  
  module procedure solve_temp_sub
    integer :: ij, ir, is, ijm
    
    if (rematrix) call this%prepare_mat_temp_sub( this%j_indx(ijmstart) , this%j_indx(ijmend) )
    
    !$omp parallel do private (ir,is,ij)
    do ijm = ijmstart, ijmend, ijmstep
      ij = this%j_indx(ijm)
      
      if ( matxsol ) then
        do concurrent ( ir=2:this%nd )
          this%rtemp(ir,ijm) = this%rtemp(ir,ijm) + this%mat%temp(ij)%multipl_fn(3*(ir-1)+1, this%sol%temp(:,ijm))
        end do
      end if
      
      do concurrent ( ir=1:this%nd )
        is = 3*(ir-1)+1
        
        this%sol%temp(is  ,ijm) = this%rtemp(ir,ijm)
        this%sol%temp(is+1,ijm) = czero
        this%sol%temp(is+2,ijm) = czero
      end do
        
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,ijm) = this%rtemp(ir,ijm)
        
      call this%mat%temp(ij)%luSolve_sub( this%sol%temp(:,ijm) )
    end do
    !$omp end parallel do
    
  end procedure solve_temp_sub
  
end submodule equations_temp