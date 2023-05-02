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

  subroutine solve_temp_sub(this)
    class(T_physicalObject), intent(inout) :: this
    integer                                :: ij, ir, ir1, ijm
    
    !$omp parallel do private (ir,ir1,ij)
    do ijm = 1, this%jms
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

end submodule Equations_temp