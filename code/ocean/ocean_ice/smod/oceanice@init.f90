submodule (oceanice) init
  implicit none; contains
  
  module procedure init_oceanice_sub
    
    !! Initialize the ocean model
    call this%init_ocean_sub()
    
    !! Initialize the equations
    call this%init_eq_temp_sub()
    call this%init_eq_torr_sub()
    call this%init_eq_mech_sub()
    
    !! Initialize the matrices
    call this%prepare_mat_temp_sub( ijstart=0 , ijend=this%jmax )
    call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
    !! Initialize the non-linear terms
    allocate( this%ntemp(this%jms,2:this%nd) ); this%ntemp = czero
    allocate( this%ntorr(this%jms,2:this%nd) ); this%ntorr = czero
    allocate( this%nsph1(this%jms,2:this%nd) ); this%nsph1 = czero
    allocate( this%nsph2(this%jms,2:this%nd) ); this%nsph2 = czero
    
    !! Set the thermal bottom boundary condition
    call this%init_temp_bbnd_sub()
    
    !! Initialize the upper layer deformation and heat flux
    call this%bnd%init_layer_up_sub()
    call this%bnd%init_flux_up_sub()
    
    !! Initialize the state
    call this%init_state_sub()
    
  end procedure init_oceanice_sub
  
end submodule init