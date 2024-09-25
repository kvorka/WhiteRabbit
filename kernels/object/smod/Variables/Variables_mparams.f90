submodule (PhysicalObject) Variables_mparams
  implicit none; contains
  
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
  
end submodule Variables_mparams