submodule (physicalobject) visc
  implicit none; contains
  
  module procedure visc_r_fn
    
    if ( this%mparams%initvisc ) then
      visc_r_fn = s4pi / c2r_fn( this%mparams%visc(1,ir) )
    else
      visc_r_fn = one
    end if
    
  end procedure visc_r_fn
  
  module procedure visc_rr_fn
    
    if ( this%mparams%initvisc ) then
      visc_rr_fn = s4pi * ( this%rad_grid%cc(ir,-1) / c2r_fn( this%mparams%visc(1,ir-1) ) + &
                          & this%rad_grid%cc(ir,+1) / c2r_fn( this%mparams%visc(1,ir  ) )   )
    else
      visc_rr_fn = one
    end if
    
  end procedure visc_rr_fn
  
end submodule visc