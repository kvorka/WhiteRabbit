submodule (PhysicalObject) Equations_temp
  implicit none
  contains

  subroutine init_eq_temp_sub(this, rhs, nl)
    class(T_physicalObject), intent(inout) :: this
    logical,                 intent(in)    :: rhs, nl
    
    call this%sol%init_stemp_sub()
    call this%mat%init_mtemp_sub()

    allocate( this%flux_up(this%jms) ) ; this%flux_up = czero

    if (rhs) then
      allocate( this%rtemp(this%nd+1,this%jms) ) ; this%rtemp = czero
    end if

    if (nl) then
      allocate( this%ntemp(this%jms,2:this%nd) ) ; this%ntemp = czero
    end if
    
  end subroutine init_eq_temp_sub

  subroutine prepare_mat_temp_sub(this, ijstart)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: ijstart
    integer                                :: ij
    
    !$omp parallel do
    do ij = ijstart, this%jmax
      call this%mat%temp(ij)%fill_sub( this%matica_temp_fn(j_in=ij, a_in=  this%cf), &
                                     & this%matica_temp_fn(j_in=ij, a_in=1-this%cf)  )
    end do
    !$omp end parallel do
    
  end subroutine prepare_mat_temp_sub

  subroutine solve_temp_sub(this, ijmstart, rematrix)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: ijmstart
    logical,                 intent(in)    :: rematrix
    integer                                :: ij, ir, is, ijm
    
    if (rematrix) call this%prepare_mat_temp_sub( this%j_indx(ijmstart) )
    
    !$omp parallel do private (ir,is,ij)
    do ijm = ijmstart, this%jms
      ij = this%j_indx(ijm)
      
      do concurrent ( ir=2:this%nd )
        this%rtemp(ir,ijm) = this%rtemp(ir,ijm) + this%mat%temp(ij)%multipl_fn(3*(ir-1)+1, this%sol%temp(:,ijm))
      end do
      
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
    
  end subroutine solve_temp_sub
  
  subroutine solve_temp_deg0_sub(this, qConv)
    class(T_physicalObject), intent(inout) :: this
    real(kind=dbl),          intent(out)   :: qConv
    integer                                :: ir, is
    complex(kind=dbl),       allocatable   :: Temp(:), Temp1(:)

    allocate( Temp(this%nd+1), Temp1(this%nd+1) ); Temp = this%sol%temp_i_fn(1)
    
    do
      Temp1 = this%sol%temp_i_fn(1)
      call this%mat%temp(0)%fill_sub( this%matica_temp_fn(j_in=0, a_in=1._dbl), this%matica_temp_fn(j_in=0, a_in=0._dbl) )
    
      ir = 1
        this%sol%temp(1,1) = cs4pi
        this%sol%temp(2,1) = czero
        this%sol%temp(3,1) = czero
  
      do ir = 2, this%nd
        is = 3*(ir-1)+1
        
        this%sol%temp(is  ,1) = Temp(ir) / this%dt + this%ntemp(1,ir) + this%htide_fn(ir,1)
        this%sol%temp(is+1,1) = czero
        this%sol%temp(is+2,1) = czero
      end do
  
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,1) = czero
  
      call this%mat%temp(0)%luSolve_sub(this%sol%temp(:,1))
  
      if ( maxval(abs(this%sol%temp_i_fn(1) - Temp1)/abs(Temp1)) < 1e-5) exit       
    end do

    deallocate( Temp, Temp1 )

    qConv = c2r_fn( -this%sol%flux_fn(1,1,1) / sqrt(4*pi) )
    
  end subroutine solve_temp_deg0_sub

end submodule Equations_temp