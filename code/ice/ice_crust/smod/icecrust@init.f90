submodule (icecrust) init
  implicit none; contains
  
  module procedure init_iceCrust_sub
    integer :: ijm
    
    call this%init_ice_sub(jmax_in = jmax_ice, rheol_in = 'viscel', n_iter = n_iter_ice)
    
    this%cf = one
    
    call this%init_eq_temp_sub( rflux=.true. )
      allocate( this%ntemp(this%jms,2:this%nd)   )
      allocate( this%nflux(3,this%jms,this%nd) )
      this%ntemp = czero
      this%nflux = czero
      
    call this%init_eq_mech_sub()
    
    call this%tdheat%init_sub(this%nd, this%jms)
    call this%bnd%init_temp_up_sub()
      call this%set_surfTemp_sub()

    open( unit=1, file='data/surf-temp.dat', status='new', action='write' )
      ijm = 1
        write(1,*) 1, this%bnd%temp_up(1) * ( this%Td - this%Tu ) + this%Tu * s4pi
        
      do ijm = 2, this%jms
        write(1,*) ijm, this%bnd%temp_up(ijm) * ( this%Td - this%Tu )
      end do
    close(1)
    
    call this%mparams%init_visc_radial_sub()
    call this%mparams%init_cp_sub()
    call this%mparams%init_lambda_sub()
    call this%mparams%init_alpha_radial_sub()
    
    call this%tides%init_sub()
    
  end procedure init_iceCrust_sub
  
end submodule init