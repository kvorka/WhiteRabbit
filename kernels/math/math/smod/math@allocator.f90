submodule (math) allocator
  implicit none; contains
  
  module procedure alloc_aligned_sub
    
    c_arr = malloc( alig, n * size_d )
    call c_f_pointer( c_arr, f_arr, [n] )
    
  end procedure alloc_aligned_sub
  
  module procedure free_aligned_sub
    
    nullify( f_arr )
    call free( c_arr )
    
  end procedure free_aligned_sub
  
end submodule allocator