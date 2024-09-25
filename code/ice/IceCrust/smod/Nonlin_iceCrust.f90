submodule (IceCrustMod) Nonlin_iceCrust
  implicit none; contains
  
  module subroutine mvgradT_sub(this, ir, mvgradT)
    class(T_iceCrust), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: mvgradT(:)
    complex(kind=dbl), allocatable :: v(:), T(:), gradT(:)
    
    allocate( v(this%jmv), T(this%jms), gradT(this%jmv) )
      
      call this%v_rr_ijml_sub(ir, v)
      call this%gradT_rr_ijml_sub(ir, T, gradT, -1)
      
      call this%lat_grid%vcvv_sub( v, gradT, mvgradT )
    
    deallocate( v, T, gradT )
    
  end subroutine mvgradT_sub
  
  module subroutine cpdivq_sub(this, ir, cpdivq)
    class(T_iceCrust), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: cpdivq(:)
    integer                        :: ijm
    complex(kind=dbl), allocatable :: divq(:), cp(:)
    
    allocate( divq(this%jms), cp(this%jms) )
      
      call this%divq_rr_ijm_sub(ir, divq, -1)
      call this%varcp_rr_ijm_sub(ir, cp)
      
      call this%lat_grid%vcss_sub( cp, divq, cpdivq )
      
    deallocate( divq )
    
  end subroutine cpdivq_sub
  
  module subroutine mvgradT_cpdivq_sub(this, ir, mvgradT)
    class(T_iceCrust), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: mvgradT(:)
    integer                        :: ijm
    complex(kind=dbl), allocatable :: v(:), T(:), gradT(:), divq(:), cp(:)
    
    allocate( v(this%jmv), T(this%jms), gradT(this%jmv), divq(this%jms), cp(this%jms) )
      
      call this%v_rr_ijml_sub(ir, v)
      call this%gradT_rr_ijml_sub(ir, T, gradT, -1)
      call this%divq_rr_ijm_sub(ir, divq, -1)
      call this%varcp_rr_ijm_sub(ir, cp)
      
      call this%lat_grid%vcss_add_vcvv_sub( cp, divq, v, gradT, mvgradT )
      
    deallocate( v, T, gradT, divq, cp )
    
  end subroutine mvgradT_cpdivq_sub
  
end submodule Nonlin_iceCrust