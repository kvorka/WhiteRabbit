submodule(PhysicalObject) TimeStepping
  implicit none ; contains
  
  module pure real(kind=dbl) function velc_crit_fn(this)
    class(T_physicalObject), intent(in) :: this
    integer                             :: ir
    real(kind=dbl)                      :: L2_rvelc, L2_zvelc
    complex(kind=dbl),      allocatable :: velc_jml_ir(:)
    
    velc_crit_fn = huge(0._dbl)

    allocate( velc_jml_ir(this%jmv) ) ; velc_jml_ir = czero
    
    do ir = 2, this%nd
      call this%sol%velocity_jml_sub(ir, velc_jml_ir)
        L2_rvelc = snorm_fn(this%jmax, ervs_fn(this%jmax, velc_jml_ir))
        L2_zvelc = vnorm_fn(this%jmax, ervv_fn(this%jmax, velc_jml_ir))
      
      velc_crit_fn = min( velc_crit_fn, this%rad_grid%r(ir) - this%rad_grid%r(ir-1) / L2_rvelc, &
                        &                          this%rad_grid%rr(ir) / this%jmax / L2_zvelc  )
    end do
    
    deallocate( velc_jml_ir )
    
  end function velc_crit_fn
  
end submodule TimeStepping