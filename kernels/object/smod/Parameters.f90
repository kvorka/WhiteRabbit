submodule(PhysicalObject) Parameters
  implicit none ; contains
  
  module pure real(kind=dbl) function lambda_fn(this, i)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    
    lambda_fn = 1._dbl
    
  end function lambda_fn
  
  module pure real(kind=dbl) function cp_fn(this, i)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    
    cp_fn = 1._dbl
    
  end function cp_fn
  
  module pure real(kind=dbl) function visc_fn(this, i)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    
    visc_fn = 1._dbl
    
  end function visc_fn
  
  module pure real(kind=dbl) function alpha_fn(this, i)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    
    alpha_fn = 1._dbl
    
  end function alpha_fn
  
end submodule Parameters