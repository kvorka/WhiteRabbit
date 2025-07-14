submodule (boundaries) inittemp
  implicit none; contains
  
  module procedure init_temp_up_sub
    
    allocate( this%temp_up(this%jms) )
      call zero_carray_sub( this%jms, this%temp_up )
    
  end procedure init_temp_up_sub
  
  module procedure init_flux_up_sub
    
    allocate( this%flux_up(this%jms) )
      call zero_carray_sub( this%jms, this%flux_up )
    
  end procedure init_flux_up_sub
  
  module procedure init_flux_dn_sub
    
    allocate( this%flux_dn(this%jms) )
      call zero_carray_sub( this%jms, this%flux_dn )
    
  end procedure init_flux_dn_sub
  
end submodule inittemp