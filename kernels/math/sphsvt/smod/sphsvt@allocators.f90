submodule (sphsvt) allocators
  implicit none; contains
  
  module procedure allocate_scalars_sub
    
    allocate( cscal(ns*this%jms2) )
      call zero_carray_sub( ns*this%jms2, cscal(1) )
    
  end procedure allocate_scalars_sub
  
  module procedure allocate_vectors_sub
    
    allocate( cvec(nv*this%jmv1) )
      call zero_carray_sub( nv*this%jmv1, cvec(1) )
    
  end procedure allocate_vectors_sub
  
end submodule allocators