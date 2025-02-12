submodule (math) allocator2
  implicit none; contains
  
  module procedure alloc_aligned2d_sub
    
    c_arr = malloc( alig, n1 * n2 * size_d )
    call c_f_pointer( c_arr, f_arr, [n1,n2] )
    
  end procedure alloc_aligned2d_sub
  
  module procedure free_aligned2d_sub
    
    nullify( f_arr )
    call free( c_arr )
    
  end procedure free_aligned2d_sub
  
end submodule allocator2