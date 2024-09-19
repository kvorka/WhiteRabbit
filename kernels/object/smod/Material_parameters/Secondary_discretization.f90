submodule(PhysicalObject) Secondary_discretization
  implicit none ; contains
  
  module pure real(kind=dbl) function lambda_rr_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    lambda_rr_fn = one
    
  end function lambda_rr_fn
  
  module pure real(kind=dbl) function cp_rr_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    cp_rr_fn = one
    
  end function cp_rr_fn
  
  module pure real(kind=dbl) function visc_rr_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    visc_rr_fn = one
    
  end function visc_rr_fn
  
  module pure real(kind=dbl) function alpha_rr_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    alpha_rr_fn = one
    
  end function alpha_rr_fn
  
end submodule Secondary_discretization