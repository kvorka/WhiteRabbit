submodule (PhysicalObject) MP_visc
  implicit none; contains
  
  module pure real(kind=dbl) function visc_r_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initvisc ) then
      visc_r_fn = s4pi / c2r_fn( this%mparams%visc(1,ir) )
    else
      visc_r_fn = one
    end if
    
  end function visc_r_fn
  
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
  
end submodule MP_visc