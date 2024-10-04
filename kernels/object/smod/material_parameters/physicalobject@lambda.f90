submodule (physicalobject) lambda
  implicit none; contains
  
  module pure real(kind=dbl) function lambda_r_fn(this, ir)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    
    if ( this%mparams%initlambda ) then
      lambda_r_fn = s4pi / c2r_fn( this%mparams%lambda(1,ir) )
    else
      lambda_r_fn = one
    end if
    
  end function lambda_r_fn
  
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
  
  pure module subroutine varlambda_r_ijm_sub(this, ir, varlambda)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: varlambda(:)
    integer                              :: ijm
    
    if ( this%mparams%initlambda .and. .not. this%mparams%lambda_radial ) then
      !ijm = 1
        varlambda(1) = czero
      
      do concurrent ( ijm = 2:this%jms )
        varlambda(ijm) = this%mparams%lambda(ijm,ir)
      end do
    
    else
      call zero_carray_sub( this%jms, varlambda )
      
    end if
    
  end subroutine varlambda_r_ijm_sub
  
end submodule lambda