submodule (lege_poly) allocators
  implicit none; contains
  
  module procedure allocate_rscalars_sub
    integer :: n
    
    n = 4*ns*this%nrma
    
    c_rscal = malloc( 32, n*size_d )
    call c_f_pointer( c_rscal, rscal, [n] )
    call zero_rarray_sub( n, rscal )
    
  end procedure allocate_rscalars_sub
  
end submodule allocators