submodule (PhysicalObject) Equations_temp
  implicit none
  contains

  subroutine init_eq_temp_sub(this, rhs, nl)
    class(T_physicalObject), intent(inout) :: this
    logical,                 intent(in)    :: rhs, nl
    integer                                :: j
    
    call this%sol%init_stemp_sub()
    call this%mat%init_mtemp_sub()

    allocate( this%flux_up(this%jms) ) ; this%flux_up = czero

    if (rhs) then
      allocate( this%rtemp(this%nd+1,this%jms) ) ; this%rtemp = czero

      this%rtemp(1,1) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
    end if

    if (nl) then
      allocate( this%ntemp(this%jms,2:this%nd) ) ; this%ntemp = czero
    end if

    do j=0, this%jmax
      call this%mat%temp(j)%fill_sub( this%matica_temp_fn(j,+0.6_dbl), this%matica_temp_fn(j,-0.4_dbl) )
    end do
    
  end subroutine init_eq_temp_sub

  subroutine solve_temp_sub(this, ijmstart)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: ijmstart
    integer                                :: ij, ir, ir1, ijm
    
    !$omp parallel do private (ir,ir1,ij)
    do ijm = ijmstart, this%jms
      ij = this%j_indx(ijm)
      
      do concurrent ( ir=2:this%nd )
        this%rtemp(ir,ijm) = this%rtemp(ir,ijm) + this%mat%temp(ij)%multipl_fn(3*(ir-1)+1, this%sol%temp(:,ijm))
      end do
      
      do concurrent ( ir=1:this%nd )
        ir1 = 3*(ir-1)+1
        
        this%sol%temp(ir1  ,ijm) = this%rtemp(ir,ijm)
        this%sol%temp(ir1+1,ijm) = czero
        this%sol%temp(ir1+2,ijm) = czero
      end do
        
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,ijm) = this%rtemp(ir,ijm)
        
      call this%mat%temp(ij)%luSolve_sub( this%sol%temp(:,ijm) )
    end do
    !$omp end parallel do
    
  end subroutine solve_temp_sub
  
  subroutine solve_temp_deg0_sub(this)
    class(T_physicalObject), intent(inout) :: this
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
    
  end subroutine solve_temp_deg0_sub

end submodule Equations_temp