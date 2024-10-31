submodule (icecrust) output
  implicit none; contains
  
  module procedure vypis_iceCrust_sub
      
    call this%vypis_sub(8, 'data/data_ice_topo' , 'topo' )
    call this%vypis_sub(8, 'data/data_ice_shape', 'shape')
    
    this%poc = this%poc + 1
    
  end procedure vypis_iceCrust_sub
  
end submodule output