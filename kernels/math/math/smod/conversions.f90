submodule (math) conversions
  implicit none; contains
  
  module elemental function int2str_fn(n) result(str)
    integer,          intent(in) :: n
    character(len=10)            :: str
    
    write(str,'(1I4)') n
    
  end function int2str_fn
  
  module elemental real(kind=dbl) function i2r_fn(ix)
    integer, intent(in) :: ix
    
    i2r_fn = real(ix, kind=dbl)
    
  end function i2r_fn
  
  module elemental complex(kind=dbl) function r2c_fn(x)
    real(kind=dbl), intent(in) :: x
    
    r2c_fn = cmplx(x, zero, kind=dbl)
    
  end function r2c_fn
  
  module elemental real(kind=dbl) function c2r_fn(cx)
    complex(kind=dbl), intent(in) :: cx
    
    c2r_fn = real(cx, kind=dbl)
    
  end function c2r_fn
  
end submodule conversions