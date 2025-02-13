submodule (lege_poly) allocators
  implicit none; contains
  
  module procedure allocate_rscalars_sub
    
    allocate( rscal(4*ns*this%nrma) )
      call zero_rarray_sub( 4*ns*this%nrma, rscal )
    
  end procedure allocate_rscalars_sub
  
end submodule allocators