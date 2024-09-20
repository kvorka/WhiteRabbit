submodule (IceTidesMod) EE_iceTides
  implicit none; contains
  
  module subroutine EE_mech_iceTides_sub(this)
    class(T_iceTides), intent(inout) :: this
    integer                          :: ijm
    
    do concurrent ( ijm = 4:6:2 )
      this%rsph1(          1,ijm) = -this%bnd%u_dn(ijm) + this%Vdelta_fn(1,ijm)
      this%rsph1(2:this%nd+1,ijm) = czero
      
      this%rsph2(1:this%nd,ijm) = czero
      this%rsph2(this%nd+1,ijm) = -this%bnd%u_up(ijm) + this%Vdelta_fn(this%nd,ijm)
    end do
    
    call this%solve_mech_sub( ijmstart=4, ijmend=6, ijmstep=2, rematrix=.true., matxsol=.false. )
    
    do concurrent ( ijm = 4:6:2 )
      this%bnd%u_dn(ijm) = this%bnd%u_dn(ijm) + this%vr_r_fn(1      ,ijm) * this%dt
      this%bnd%u_up(ijm) = this%bnd%u_up(ijm) + this%vr_r_fn(this%nd,ijm) * this%dt
    end do
    
  end subroutine EE_mech_iceTides_sub
  
end submodule EE_iceTides