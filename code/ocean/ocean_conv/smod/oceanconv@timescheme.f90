submodule (oceanconv) timescheme
  implicit none; contains
  
  module procedure time_scheme_oceanConv_sub
    integer :: ir, ijm
    
    !ir = 1 ; ijm = 1
      this%rtemp(1,1) = cs4pi
    
    !$omp parallel
    !$omp do
    do ijm = 2, this%jms
      this%rsph2(this%nd+1,ijm) = -this%St * this%qr_r_fn(this%nd,ijm)
    end do
    !$omp end do

    !$omp do collapse (2) schedule (dynamic, ocean_chunk_size)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtemp(ir,ijm) = (1-this%ab) * this%ntemp(ijm,ir)
        this%rtorr(ir,ijm) = (1-this%ab) * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = (1-this%ab) * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = (1-this%ab) * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end do
    !$omp end parallel
    
    call this%fullnl_sub()
    
    !$omp parallel do collapse (2) schedule (dynamic, ocean_chunk_size)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtemp(ir,ijm) = this%rtemp(ir,ijm) + this%ab * this%ntemp(ijm,ir)
        this%rtorr(ir,ijm) = this%rtorr(ir,ijm) + this%ab * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + this%ab * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + this%ab * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end parallel do
    
    call this%solve_temp_sub( ijmstart=1 , ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.true. )
    call this%solve_torr_sub( ijmstart=2 , ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.true. )
    call this%solve_mech_sub( ijmstart=2 , ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.true. )
    
    if (this%mechanic_bnd == 'frees') call this%global_rotation_sub()
    
  end procedure time_scheme_oceanConv_sub
  
end submodule timescheme