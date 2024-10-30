submodule (solution) init
  implicit none; contains
  
  module procedure init_solution_sub
    
    this%nd   = nd
    this%jmax = jmax
    this%jms  =     jmax * (jmax+1) / 2 + jmax   + 1
    this%jmv  = 3*( jmax * (jmax+1) / 2 + jmax ) + 1
    this%jmt  = 5*( jmax * (jmax+1) / 2 + jmax ) + 1
    
    this%inittemp = .false.
    this%initsfer = .false.
    this%inittorr = .false.
    
  end procedure init_solution_sub
  
  module procedure deallocate_solution_sub
    
    if ( allocated(this%temp) ) deallocate( this%temp )
    if ( allocated(this%torr) ) deallocate( this%torr )
    if ( allocated(this%mech) ) deallocate( this%mech )
    
  end procedure deallocate_solution_sub
  
  module procedure init_stemp_sub
    
    allocate( this%temp(3*this%nd+1, this%jms) )
    
    this%inittemp = .true.
    this%temp     = czero
    
  end procedure init_stemp_sub
  
  module procedure init_storr_sub
      
    allocate( this%torr(3*this%nd+1, this%jms) )
    
    this%inittorr = .true.
    this%torr = czero
    
  end procedure init_storr_sub
  
  module procedure init_smech_sub
    
    allocate( this%mech(6*this%nd+2,this%jms) )
    
    this%initsfer = .true.
    this%mech     = czero
    
  end procedure init_smech_sub
  
end submodule init