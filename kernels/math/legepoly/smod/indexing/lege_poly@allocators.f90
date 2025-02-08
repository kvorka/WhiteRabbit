submodule (lege_poly) allocators
  implicit none; contains
  
  module procedure allocate_cscalars_sub
    
    allocate( cscal(2*ns*this%nrma) ); cscal = czero
    
  end procedure allocate_cscalars_sub
  
end submodule allocators