submodule(PhysicalObject) Main_discretization
  implicit none ; contains
  
  module pure real(kind=dbl) function lambda_r_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    lambda_r_fn = one
    
  end function lambda_r_fn
  
  module pure real(kind=dbl) function cp_r_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    cp_r_fn = one
    
  end function cp_r_fn
  
  module pure real(kind=dbl) function visc_r_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    visc_r_fn = one
    
  end function visc_r_fn
  
  module pure real(kind=dbl) function alpha_r_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    alpha_r_fn = one
    
  end function alpha_r_fn
  
end submodule Main_discretization