submodule (physicalobject) cp
  implicit none; contains
  
  module procedure cp1_r_fn
    
    if ( this%mparams%initcp ) then
      cp1_r_fn = ( this%rad_grid%c(ir,-1) * c2r_fn( this%mparams%cp(1,ir  ) ) + &
                 & this%rad_grid%c(ir,+1) * c2r_fn( this%mparams%cp(1,ir+1) )   ) / s4pi
    else
      cp1_r_fn = one
    end if
    
  end procedure cp1_r_fn
  
  module procedure cp1_rr_fn
    
    if ( this%mparams%initcp ) then
      cp1_rr_fn = c2r_fn( this%mparams%cp(1,ir) ) / s4pi
    else
      cp1_rr_fn = one
    end if
    
  end procedure cp1_rr_fn
  
  module procedure varcp1_rr_ijm_sub
    integer :: ijm
    
    if ( this%mparams%initcp .and. .not. this%mparams%cp_radial ) then
      !ijm = 1
        varcp(1) = czero
      
      do concurrent ( ijm = 2:this%jms )
        varcp(ijm) = this%mparams%cp(ijm,ir)
      end do
    
    else
      call zero_carray_sub( this%jms, varcp )
      
    end if
    
  end procedure varcp1_rr_ijm_sub
  
end submodule cp