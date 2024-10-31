submodule (icecrust) nonlin
  implicit none; contains
  
  module procedure mvgradT_sub
    complex(kind=dbl), allocatable :: v(:), T(:), gradT(:)
    
    allocate( v(this%jmv), T(this%jms), gradT(this%jmv) )
      
      call this%v_rr_ijml_sub(ir, v)
      call this%gradT_rr_ijml_sub(ir, T, gradT, -1)
      
      call this%lat_grid%vcvv_sub( v, gradT, mvgradT )
    
    deallocate( v, T, gradT )
    
  end procedure mvgradT_sub
  
  module procedure cpdivq_sub
    integer                        :: ijm
    complex(kind=dbl), allocatable :: divq(:), cp(:)
    
    allocate( divq(this%jms), cp(this%jms) )
      
      call this%divq_rr_ijm_sub(ir, divq, -1)
      call this%varcp_rr_ijm_sub(ir, cp)
      
      call this%lat_grid%vcss_sub( cp, divq, cpdivq )
      
    deallocate( divq )
    
  end procedure cpdivq_sub
  
  module procedure mvgradT_cpdivq_sub
    integer                        :: ir
    complex(kind=dbl), allocatable :: v(:), T(:), gradT(:), divq(:), cp(:)
    
    allocate( v(this%jmv), T(this%jms), gradT(this%jmv), divq(this%jms), cp(this%jms) )
    
    !$omp parallel do private (v, T, gradT, divq, cp)
    do ir = 2, this%nd
      call this%v_rr_ijml_sub(ir, v)
      call this%gradT_rr_ijml_sub(ir, T, gradT, -1)
      call this%divq_rr_ijm_sub(ir, divq, -1)
      call this%varcp_rr_ijm_sub(ir, cp)
      
      call this%lat_grid%vcss_add_vcvv_sub( cp, divq, v, gradT, this%ntemp(:,ir) )
    end do
    !$omp end parallel do
    
    deallocate( v, T, gradT, divq, cp )
    
  end procedure mvgradT_cpdivq_sub
  
end submodule nonlin