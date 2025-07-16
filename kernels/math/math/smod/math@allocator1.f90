submodule (math) allocator1
  implicit none; contains
  
  module procedure alloc_aligned1d_sub
    
    c_arr = malloc( alig, n * int( c_sizeof(0._dbl) ) )
    call c_f_pointer( c_arr, f_arr, [n] )
    
  end procedure alloc_aligned1d_sub
  
  module procedure free_aligned1d_sub
    
    nullify( f_arr )
    call free( c_arr )
    
  end procedure free_aligned1d_sub
  
end submodule allocator1