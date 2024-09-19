submodule(PhysicalObject) Variables_tideheat
  implicit none ; contains
  
  module pure complex(kind=dbl) function htide_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    
    htide_fn = this%Ds/this%Ra * ( this%rad_grid%cc(ir,-1) * this%htide(ir-1,ijm) + &
                                 & this%rad_grid%cc(ir,+1) * this%htide(ir  ,ijm)   ) / this%cp_rr_fn(ir)
    
  end function htide_fn
  
end submodule Variables_tideheat