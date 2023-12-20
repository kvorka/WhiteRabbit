submodule(PhysicalObject) TimeStepping
  implicit none ; contains
  
  module real(kind=dbl) function velc_crit_fn(this)
    class(T_physicalObject), intent(in) :: this
    integer                             :: ir
    real(kind=dbl)                      :: L2_rvelc, L2_zvelc
    real(kind=dbl),         allocatable :: crit(:)
    complex(kind=dbl),      allocatable :: velc_jml_ir(:)
    
    allocate( velc_jml_ir(this%jmv) , crit(2:this%nd) )
    
    !$omp parallel do private( velc_jml_ir, L2_rvelc, L2_zvelc )
    do ir = 2, this%nd
      call this%sol%velocity_jml_sub(ir, velc_jml_ir)
        L2_rvelc = max( snorm_fn(this%jmax, ervs_fn(this%jmax, velc_jml_ir)) , 1d-15 )
        L2_zvelc = max( vnorm_fn(this%jmax, ervv_fn(this%jmax, velc_jml_ir)) , 1d-15 )
      
      crit(ir) = min( ( this%rad_grid%r(ir) - this%rad_grid%r(ir-1) ) / L2_rvelc, this%rad_grid%rr(ir) / this%jmax / L2_zvelc )
    end do
    !$omp end parallel do
    
    velc_crit_fn = minval( crit )
    
    deallocate( velc_jml_ir , crit )
    
  end function velc_crit_fn
  
end submodule TimeStepping