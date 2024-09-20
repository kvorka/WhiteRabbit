submodule(PhysicalObject) Variables_tideheat
  implicit none ; contains
  
  module pure complex(kind=dbl) function htide_r_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    
    htide_r_fn = this%Ds/this%Ra * this%htide(ir,ijm) / this%cp_r_fn(ir)
    
  end function htide_r_fn
  
  module pure subroutine htide_ir_ijm_sub(this, htide)
    class(T_physicalObject), intent(in)    :: this
    complex(kind=dbl),       intent(inout) :: htide(:,:)
    integer                                :: ijm, ir
    
    do concurrent ( ijm = 1:size(this%htide, dim=2), ir = 1:this%nd )
      htide(ir,ijm) = this%htide(ir,ijm)
    end do
    
  end subroutine htide_ir_ijm_sub
  
  module pure complex(kind=dbl) function htide_rr_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    
    htide_rr_fn = this%Ds/this%Ra * ( this%rad_grid%cc(ir,-1) * this%htide(ir-1,ijm) + &
                                    & this%rad_grid%cc(ir,+1) * this%htide(ir  ,ijm)   ) / this%cp_rr_fn(ir)
    
  end function htide_rr_fn
  
end submodule Variables_tideheat