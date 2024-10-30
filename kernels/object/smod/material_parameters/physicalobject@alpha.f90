submodule (physicalobject) alpha
  implicit none; contains
  
  module procedure alpha_r_fn
    
    if ( this%mparams%initalpha ) then
      alpha_r_fn = ( this%rad_grid%c(ir,-1) * c2r_fn( this%mparams%alpha(1,ir  ) ) + &
                   & this%rad_grid%c(ir,+1) * c2r_fn( this%mparams%alpha(1,ir+1) )   ) / s4pi
    else
      alpha_r_fn = one
    end if
    
  end procedure alpha_r_fn
  
  module procedure alpha_rr_fn
    
    if ( this%mparams%initalpha ) then
      alpha_rr_fn = c2r_fn( this%mparams%alpha(1,ir) ) / s4pi
    else
      alpha_rr_fn = one
    end if
    
  end procedure alpha_rr_fn
  
end submodule alpha