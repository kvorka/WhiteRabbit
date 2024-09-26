submodule (PhysicalObject) MP_cp
  implicit none; contains
  
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
  
  module pure real(kind=dbl) function cp_rr_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initcp ) then
      cp_rr_fn = s4pi / c2r_fn( this%mparams%cp(1,ir) )
    else
      cp_rr_fn = one
    end if
    
  end function cp_rr_fn
  
  pure module subroutine varcp_rr_ijm_sub(this, ir, varcp)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: varcp(:)
    integer                              :: ijm
    
    if ( this%mparams%initcp .and. .not. this%mparams%cp_radial ) then
      !ijm = 1
        varcp(1) = czero
      
      do concurrent ( ijm = 2:this%jms )
        varcp(ijm) = this%mparams%cp(ijm,ir)
      end do
    
    else
      call zero_carray_sub( this%jms, varcp )
      
    end if
    
  end subroutine varcp_rr_ijm_sub
  
end submodule MP_cp