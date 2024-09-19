submodule(PhysicalObject) Secondary_discretization
  implicit none ; contains
  
  module pure real(kind=dbl) function lambda_rr_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initlambda ) then
      lambda_rr_fn = s4pi * ( this%rad_grid%cc(ir,-1) / c2r_fn( this%mparams%lambda(1,ir-1) ) + &
                            & this%rad_grid%cc(ir,+1) / c2r_fn( this%mparams%lambda(1,ir  ) )   )
    else
      lambda_rr_fn = one
    end if
    
  end function lambda_rr_fn
  
  module pure real(kind=dbl) function cp_rr_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initcp ) then
      cp_rr_fn = s4pi / c2r_fn( this%mparams%cp(1,ir) )
    else
      cp_rr_fn = one
    end if
    
  end function cp_rr_fn
  
  module pure real(kind=dbl) function visc_rr_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initvisc ) then
      visc_rr_fn = s4pi * ( this%rad_grid%cc(ir,-1) / c2r_fn( this%mparams%visc(1,ir-1) ) + &
                          & this%rad_grid%cc(ir,+1) / c2r_fn( this%mparams%visc(1,ir  ) )   )
    else
      visc_rr_fn = one
    end if
    
  end function visc_rr_fn
  
  module pure real(kind=dbl) function alpha_rr_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initalpha ) then
      alpha_rr_fn = c2r_fn( this%mparams%alpha(1,ir) ) / s4pi
    else
      alpha_rr_fn = one
    end if
    
  end function alpha_rr_fn
  
end submodule Secondary_discretization