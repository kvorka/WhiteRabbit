submodule (PhysicalObject) Equations_mech
  implicit none
  contains

  subroutine init_eq_mech_sub(this, rhs, nl)
    class(T_physicalObject), intent(inout) :: this
    logical,                 intent(in)    :: rhs, nl
    integer                                :: j
    
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

    do j=1, this%jmax
      call this%mat%mech(j)%fill_sub( this%matica_mech_fn(j,+0.6_dbl), this%matica_mech_fn(j,-0.4_dbl) )
    end do

  end subroutine init_eq_mech_sub

  subroutine solve_mech_sub(this)
    class(T_physicalObject), intent(inout) :: this
    integer                                :: ij, ijm, ir, ir1

    !$omp parallel do private (ir,ir1,ij)
    do ijm = 2, this%jms
      ij = this%j_indx(ijm)

      do concurrent ( ir=2:this%nd )
        this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + this%mat%mech(ij)%multipl_fn(6*(ir-1)+1,this%sol%mech(:,ijm))
        this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + this%mat%mech(ij)%multipl_fn(6*(ir-1)+2,this%sol%mech(:,ijm))
      end do
      
      do concurrent ( ir=1:this%nd )
        ir1 = 6*(ir-1)+1
        
        this%sol%mech(ir1  ,ijm) = this%rsph1(ir,ijm)
        this%sol%mech(ir1+1,ijm) = this%rsph2(ir,ijm)
        this%sol%mech(ir1+2,ijm) = czero
        this%sol%mech(ir1+3,ijm) = czero
        this%sol%mech(ir1+4,ijm) = czero
        this%sol%mech(ir1+5,ijm) = czero
      end do
        
      ir = this%nd+1
        this%sol%mech(6*this%nd+1,ijm) = czero
        this%sol%mech(6*this%nd+2,ijm) = czero
        
      call this%mat%mech(ij)%luSolve_sub( this%sol%mech(:,ijm) )
    end do
    !$omp end parallel do

  end subroutine solve_mech_sub

end submodule Equations_mech