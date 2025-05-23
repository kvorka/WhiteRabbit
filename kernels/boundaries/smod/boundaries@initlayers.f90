submodule (boundaries) initlayers
  implicit none; contains
  
  module procedure init_layers_sub
    
    allocate( this%u_up(this%jms) ); this%u_up = czero
    allocate( this%u_dn(this%jms) ); this%u_dn = czero
    allocate( this%u_I2(this%jms) ); this%u_I2 = czero
    allocate( this%u_C(this%jms)  ); this%u_C  = czero
    
    allocate( this%t_dn(this%jms) ); this%t_dn = czero
    allocate( this%t_up(this%jms) ); this%t_up = czero
    
    allocate( this%v_up(this%jms) ); this%v_up = czero
    allocate( this%v_dn(this%jms) ); this%v_dn = czero
    
  end procedure init_layers_sub
  
  module procedure init_layer_up_sub
    
    allocate( this%t_up(this%jms) ); this%t_up = czero
    allocate( this%u_up(this%jms) ); this%u_up = czero
    
  end procedure init_layer_up_sub
  
end submodule initlayers