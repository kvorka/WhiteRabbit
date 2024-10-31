submodule (oceantides) timescheme
  implicit none; contains
  
  module procedure time_scheme_oceanTides_sub
    integer :: ir, ijm
    
    !$omp parallel do collapse (2)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtorr(ir,ijm) = (1-this%ab) * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = (1-this%ab) * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = (1-this%ab) * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end parallel do
    
    if ( this%noharm ) then
      call this%coriolis_sub()
    else
      call this%coriolis_vgradv_sub()
    end if
    
    !$omp parallel do collapse (2)
    do ijm = 1, this%jms
      do ir = 2, this%nd
        this%rtorr(ir,ijm) = this%rtorr(ir,ijm) + this%ab * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + this%ab * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + this%ab * this%nsph2(ijm,ir)
      end do
    end do
    !$omp end parallel do
    
    !$omp parallel do
    do ijm = 2, this%jms
      this%rtorr(1,ijm) = czero
      this%rsph1(1,ijm) = czero
      this%rsph2(1,ijm) = czero
      
      this%rtorr(this%nd+1,ijm) = czero
      
      if (ijm == 4) then
        this%rsph1(this%nd+1,ijm) = this%v201(this%k_of_period)
        this%rsph2(this%nd+1,ijm) = this%v203(this%k_of_period)
      else if (ijm == 6) then
        this%rsph1(this%nd+1,ijm) = this%v221(this%k_of_period)
        this%rsph2(this%nd+1,ijm) = this%v223(this%k_of_period)
      else
        this%rsph1(this%nd+1,ijm) = czero
        this%rsph2(this%nd+1,ijm) = czero
      end if
    end do
    !$omp end parallel do
    
    call this%solve_torr_sub( ijmstart=2 , ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.true. )
    call this%solve_mech_sub( ijmstart=2 , ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.true. )
    
  end procedure time_scheme_oceanTides_sub
  
end submodule timescheme