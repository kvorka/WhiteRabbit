submodule (oceantides) init
  implicit none; contains
  
  module procedure init_oceanTides_sub
    
    call this%init_ocean_sub()
    
    this%dt                = 2 * pi / this%n_iter
    this%number_of_periods = 0
    
    call this%init_eq_torr_sub()
    call this%init_eq_mech_sub()
    
    allocate( this%ntorr(this%jms,2:this%nd) ); this%ntorr = czero
    allocate( this%nsph1(this%jms,2:this%nd) ); this%nsph1 = czero
    allocate( this%nsph2(this%jms,2:this%nd) ); this%nsph2 = czero
    
    call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
    call this%init_ubnd_sub()
    
  end procedure init_oceanTides_sub
  
end submodule init