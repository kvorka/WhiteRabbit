submodule(PhysicalObject) Main_discretization
  implicit none ; contains
  
  module pure real(kind=dbl) function lambda_r_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initlambda ) then
      lambda_r_fn = s4pi / c2r_fn( this%mparams%lambda(1,ir) )
    else
      lambda_r_fn = one
    end if
    
  end function lambda_r_fn
  
  module pure real(kind=dbl) function cp_r_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initcp ) then
      cp_r_fn = s4pi * ( this%rad_grid%c(ir,-1) / c2r_fn( this%mparams%cp(1,ir  ) ) + &
                       & this%rad_grid%c(ir,+1) / c2r_fn( this%mparams%cp(1,ir+1) )   )
    else
      cp_r_fn = one
    end if
    
  end function cp_r_fn
  
  module pure real(kind=dbl) function visc_r_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initvisc ) then
      visc_r_fn = s4pi / c2r_fn( this%mparams%visc(1,ir) )
    else
      visc_r_fn = one
    end if
    
  end function visc_r_fn
  
  module pure real(kind=dbl) function alpha_r_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initalpha ) then
      alpha_r_fn = ( this%rad_grid%c(ir,-1) * c2r_fn( this%mparams%alpha(1,ir  ) ) + &
                   & this%rad_grid%c(ir,+1) * c2r_fn( this%mparams%alpha(1,ir+1) )   ) / s4pi
    else
      alpha_r_fn = one
    end if
    
  end function alpha_r_fn
  
end submodule Main_discretization