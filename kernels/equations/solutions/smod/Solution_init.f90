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
  
  module pure subroutine deallocate_solution_sub(this)
    class(T_solution), intent(inout) :: this
    
    if ( allocated(this%temp) ) deallocate( this%temp )
    if ( allocated(this%torr) ) deallocate( this%torr )
    if ( allocated(this%mech) ) deallocate( this%mech )
    
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
  
end submodule Solution_init