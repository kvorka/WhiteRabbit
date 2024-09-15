submodule(Solution) Solution_init
  implicit none; contains
  
  module pure subroutine init_solution_sub(this, nd, jmax)
    class(T_solution), intent(inout) :: this
    integer,           intent(in)    :: nd, jmax
    
    this%nd   = nd
    this%jmax = jmax
    this%jms  =     jmax * (jmax+1) / 2 + jmax   + 1
    this%jmv  = 3*( jmax * (jmax+1) / 2 + jmax ) + 1
    this%jmt  = 5*( jmax * (jmax+1) / 2 + jmax ) + 1
    
    this%inittemp = .false.
    this%initsfer = .false.
    this%inittorr = .false.
    
  end subroutine init_solution_sub
  
  module pure subroutine nulify_solution_sub(this)
    class(T_solution), intent(inout) :: this
    
    if ( allocated(this%temp) ) this%temp = czero
    if ( allocated(this%torr) ) this%torr = czero
    if ( allocated(this%mech) ) this%mech = czero
    
    if ( allocated(this%u_dn) ) this%u_dn = czero
    if ( allocated(this%u_up) ) this%u_up = czero
    if ( allocated(this%u_I2) ) this%u_I2 = czero
    if ( allocated(this%u_C ) ) this%u_C  = czero
    
    if ( allocated(this%t_dn) ) this%t_dn = czero
    if ( allocated(this%t_up) ) this%t_up = czero
    
    if ( allocated(this%v_up) ) this%v_up = czero
    if ( allocated(this%v_dn) ) this%v_dn = czero
    
  end subroutine nulify_solution_sub
  
  module pure subroutine deallocate_solution_sub(this)
    class(T_solution), intent(inout) :: this
    
    if ( allocated(this%temp) ) deallocate( this%temp )
    if ( allocated(this%torr) ) deallocate( this%torr )
    if ( allocated(this%mech) ) deallocate( this%mech )
    
    if ( allocated(this%u_dn) ) deallocate( this%u_dn )
    if ( allocated(this%u_up) ) deallocate( this%u_up )
    if ( allocated(this%u_I2) ) deallocate( this%u_I2 )
    if ( allocated(this%u_C ) ) deallocate( this%u_C  )
    
    if ( allocated(this%t_dn) ) deallocate( this%t_dn )
    if ( allocated(this%t_up) ) deallocate( this%t_up )
    
    if ( allocated(this%v_up) ) deallocate( this%v_up )
    if ( allocated(this%v_dn) ) deallocate( this%v_dn )
    
  end subroutine deallocate_solution_sub
  
  module pure subroutine init_stemp_sub(this)
    class(T_solution), intent(inout) :: this
    
    allocate( this%temp(3*this%nd+1, this%jms) )
    
    this%inittemp = .true.
    this%temp     = czero
    
  end subroutine init_stemp_sub
  
  module pure subroutine init_storr_sub(this)
    class(T_solution), intent(inout) :: this
      
    allocate( this%torr(3*this%nd+1, this%jms) )
    
    this%inittorr = .true.
    this%torr = czero
    
  end subroutine init_storr_sub
  
  module pure subroutine init_smech_sub(this)
    class(T_solution), intent(inout) :: this
    
    allocate( this%mech(6*this%nd+2,this%jms) )
    
    this%initsfer = .true.
    this%mech     = czero
    
  end subroutine init_smech_sub
  
  module pure subroutine init_layers_sub(this)
    class(T_solution), intent(inout) :: this
    
    allocate( this%u_up(this%jms) ); this%u_up = czero
    allocate( this%u_dn(this%jms) ); this%u_dn = czero
    allocate( this%u_I2(this%jms) ); this%u_I2 = czero
    allocate( this%u_C(this%jms)  ); this%u_C  = czero
    
    allocate( this%t_dn(this%jms) ); this%t_dn = czero
    allocate( this%t_up(this%jms) ); this%t_up = czero
    
    allocate( this%v_up(this%jms) ); this%v_up = czero
    allocate( this%v_dn(this%jms) ); this%v_dn = czero
    
  end subroutine init_layers_sub
  
  module pure subroutine init_layer_u_sub(this)
    class(T_solution), intent(inout) :: this
    
    allocate( this%t_up(this%jms) ); this%t_up = czero
    allocate( this%u_up(this%jms) ); this%u_up = czero
    
  end subroutine init_layer_u_sub
  
end submodule Solution_init