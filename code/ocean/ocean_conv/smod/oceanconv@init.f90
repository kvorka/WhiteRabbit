submodule (oceanconv) init
  implicit none; contains
  
  module subroutine init_oceanConv_sub(this)
    class(T_oceanConv), intent(inout) :: this
    
    call this%init_ocean_sub()
    
    call this%init_eq_temp_sub()
      allocate( this%ntemp(this%jms,2:this%nd) )
      this%ntemp = czero
      
      call this%prepare_mat_temp_sub( ijstart=0 , ijend=this%jmax )
    
    call this%init_eq_torr_sub()
      allocate( this%ntorr(this%jms,2:this%nd) )
      this%ntorr = czero
      
      call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    
    call this%init_eq_mech_sub()
      allocate( this%nsph1(this%jms,2:this%nd), this%nsph2(this%jms,2:this%nd) )
      this%nsph1 = czero
      this%nsph2 = czero
      
      call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
    call this%init_state_sub()
    
  end subroutine init_oceanConv_sub
  
end submodule init