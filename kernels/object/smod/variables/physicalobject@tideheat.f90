submodule (physicalobject) tideheat
  implicit none ; contains
  
  module procedure htide_r_fn
    
    htide_r_fn = this%Ds/this%Ra * this%tdheat%htide(ir,ijm) / this%cp_r_fn(ir)
    
  end procedure htide_r_fn
  
  module procedure htide_ir_ijm_sub
    integer :: ijm, ir
    
    do concurrent ( ijm = 1:this%tdheat%jms, ir = 1:this%nd )
      htide(ir,ijm) = this%tdheat%htide(ir,ijm)
    end do
    
  end procedure htide_ir_ijm_sub
  
  module procedure htide_rr_fn
    
    htide_rr_fn = this%Ds/this%Ra * ( this%rad_grid%cc(ir,-1) * this%tdheat%htide(ir-1,ijm) + &
                                    & this%rad_grid%cc(ir,+1) * this%tdheat%htide(ir  ,ijm)   ) / this%cp_rr_fn(ir)
    
  end procedure htide_rr_fn
  
end submodule tideheat