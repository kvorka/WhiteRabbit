submodule (PhysicalObject) Equations_torr
  implicit none
  contains

  subroutine init_eq_torr_sub(this, rhs, nl)
    class(T_physicalObject), intent(inout) :: this
    logical,                 intent(in)    :: rhs, nl
    
    call this%sol%init_storr_sub()
    call this%mat%init_mtorr_sub()

    if (rhs) then
      allocate( this%rtorr(this%nd+1,this%jms) ) ; this%rtorr = czero
    end if

    if (nl) then
      allocate( this%ntorr(this%jms,2:this%nd) ) ; this%ntorr = czero
    end if
    
  end subroutine init_eq_torr_sub
  
  subroutine prepare_mat_torr_sub(this, ijstart)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: ijstart
    integer                                :: ij
    
    !$omp parallel do
    do ij = ijstart, this%jmax
      call this%mat%torr(ij)%fill_sub( this%matica_torr_fn(j_in=ij, a_in=this%cf  ), &
                                     & this%matica_torr_fn(j_in=ij, a_in=this%cf-1)  )
    end do
    !$omp end parallel do
    
  end subroutine prepare_mat_torr_sub

  subroutine solve_torr_sub(this, ijmstart, rematrix)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: ijmstart
    logical,                 intent(in)    :: rematrix
    integer                                :: ij, ijm, ir, ir1
    
    if (rematrix) call this%prepare_mat_torr_sub( this%j_indx(ijmstart) )
    
    !$omp parallel do private (ir,ir1,ij)
    do ijm = 2, this%jms
      ij = this%j_indx(ijm)

      do concurrent ( ir=2:this%nd )
        this%rtorr(ir,ijm) = this%rtorr(ir,ijm) + this%mat%torr(ij)%multipl_fn(3*(ir-1)+1,this%sol%torr(:,ijm))
      end do

      do concurrent ( ir=1:this%nd )
        ir1 = 3*(ir-1) + 1
        
        this%sol%torr(ir1  ,ijm) = this%rtorr(ir,ijm)
        this%sol%torr(ir1+1,ijm) = czero
        this%sol%torr(ir1+2,ijm) = czero
      end do
  
      ir = this%nd+1
        this%sol%torr(3*this%nd+1,ijm) = this%rtorr(ir,ijm)
        
      call this%mat%torr(ij)%luSolve_sub( this%sol%torr(:,ijm) )
    end do
    !$omp end parallel do

  end subroutine solve_torr_sub

end submodule Equations_torr