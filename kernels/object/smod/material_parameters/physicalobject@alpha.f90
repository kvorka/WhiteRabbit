submodule (physicalobject) alpha
  implicit none; contains
  
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
  
  module pure real(kind=dbl) function alpha_rr_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initalpha ) then
      alpha_rr_fn = c2r_fn( this%mparams%alpha(1,ir) ) / s4pi
    else
      alpha_rr_fn = one
    end if
    
  end function alpha_rr_fn
  
end submodule alpha