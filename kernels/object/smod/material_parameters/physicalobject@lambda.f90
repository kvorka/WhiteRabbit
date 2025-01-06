submodule (physicalobject) lambda
  implicit none; contains
  
  module procedure lambda_r_fn
    
    if ( this%mparams%initlambda ) then
      lambda_r_fn = c2r_fn( this%mparams%lambda(1,ir) ) / s4pi
    else
      lambda_r_fn = one
    end if
    
  end procedure lambda_r_fn
  
  module procedure lambda_rr_fn
    
    if ( this%mparams%initlambda ) then
      lambda_rr_fn = ( this%rad_grid%cc(ir,-1) * c2r_fn( this%mparams%lambda(1,ir-1) ) + &
                     & this%rad_grid%cc(ir,+1) * c2r_fn( this%mparams%lambda(1,ir  ) )   ) / s4pi
    else
      lambda_rr_fn = one
    end if
    
  end procedure lambda_rr_fn
  
  module procedure varlambda_r_ijm_sub
    integer :: ijm
    
    if ( this%mparams%initlambda .and. .not. this%mparams%lambda_radial ) then
      !ijm = 1
        varlambda(1) = czero
      
      do concurrent ( ijm = 2:this%jms )
        varlambda(ijm) = this%mparams%lambda(ijm,ir)
      end do
    
    else
      call zero_carray_sub( this%jms, varlambda )
      
    end if
    
  end procedure varlambda_r_ijm_sub
  
end submodule lambda