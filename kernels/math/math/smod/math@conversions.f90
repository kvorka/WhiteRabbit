submodule (math) conversions
  implicit none; contains
  
  module procedure int2str_fn
    
    write(str,'(1I4)') n
    
  end procedure int2str_fn
  
  module procedure i2r_fn
    
    i2r_fn = real(ix, kind=dbl)
    
  end procedure i2r_fn
  
  module procedure r2c_fn
    
    r2c_fn = cmplx(x, zero, kind=dbl)
    
  end procedure r2c_fn
  
  module procedure c2r_fn
    
    c2r_fn = real(cx, kind=dbl)
    
  end procedure c2r_fn
  
end submodule conversions