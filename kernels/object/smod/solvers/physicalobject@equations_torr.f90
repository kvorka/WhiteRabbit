submodule (physicalobject) equations_torr
  implicit none ; contains
  
  module subroutine init_eq_torr_sub(this)
    class(T_physicalObject), intent(inout) :: this
    
    call this%sol%init_storr_sub()
    call this%mat%init_mtorr_sub()
    
    allocate( this%rtorr(this%nd+1,this%jms) )
      this%rtorr = czero
    
  end subroutine init_eq_torr_sub
  
  module subroutine prepare_mat_torr_sub(this, ijstart, ijend)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: ijstart, ijend
    integer                                :: ij
    
    !$omp parallel do
    do ij = ijstart, ijend
      call this%mat%torr(ij)%fill_sub( this%mat_torr_fn(j_in=ij, a_in=this%cf  ), &
                                     & this%mat_torr_fn(j_in=ij, a_in=this%cf-1)  )
    end do
    !$omp end parallel do
    
  end subroutine prepare_mat_torr_sub
  
  module subroutine solve_torr_sub(this, ijmstart, ijmend, ijmstep, rematrix, matxsol)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: ijmstart, ijmend, ijmstep
    logical,                 intent(in)    :: rematrix, matxsol
    integer                                :: ij, ijm, ir, is
    
    if (rematrix) call this%prepare_mat_torr_sub( this%j_indx(ijmstart) , this%j_indx(ijmend) )
    
    !$omp parallel do private (ir,is,ij)
    do ijm = ijmstart, ijmend, ijmstep
      ij = this%j_indx(ijm)

      if ( matxsol ) then
        do concurrent ( ir=2:this%nd )
          this%rtorr(ir,ijm) = this%rtorr(ir,ijm) + this%mat%torr(ij)%multipl_fn(3*(ir-1)+1,this%sol%torr(:,ijm))
        end do
      end if

      do concurrent ( ir=1:this%nd )
        is = 3*(ir-1) + 1
        
        this%sol%torr(is  ,ijm) = this%rtorr(ir,ijm)
        this%sol%torr(is+1,ijm) = czero
        this%sol%torr(is+2,ijm) = czero
      end do
  
      ir = this%nd+1
        this%sol%torr(3*this%nd+1,ijm) = this%rtorr(ir,ijm)
        
      call this%mat%torr(ij)%luSolve_sub( this%sol%torr(:,ijm) )
    end do
    !$omp end parallel do

  end subroutine solve_torr_sub
  
end submodule equations_torr