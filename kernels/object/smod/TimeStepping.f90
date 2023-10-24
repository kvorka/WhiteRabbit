submodule(PhysicalObject) TimeStepping
  implicit none
  contains
  
  real(kind=dbl) pure function velc_crit_fn(this)
    class(T_physicalObject), intent(in) :: this
    integer                             :: ir
    real(kind=dbl)                      :: dr, r, L2_rvelc, L2_zvelc
    complex(kind=dbl),      allocatable :: velc_jml_ir(:), rvelc_jm_ir(:)
    
    velc_crit_fn = huge(0._dbl)

    allocate( velc_jml_ir(this%jmv) ) ; velc_jml_ir = czero
    allocate( rvelc_jm_ir(this%jms) ) ; rvelc_jm_ir = czero
    
    do ir = 2, this%nd
      dr = this%rad_grid%r(ir) - this%rad_grid%r(ir-1)
      r  = this%rad_grid%rr(ir)
      
      call this%sol%velocity_jml_sub(ir, velc_jml_ir)
        rvelc_jm_ir =               ervs_fn(this%jmax, velc_jml_ir) ; L2_rvelc = snorm_fn(this%jmax, rvelc_jm_ir)
        velc_jml_ir = velc_jml_ir - ersv_fn(this%jmax, rvelc_jm_ir) ; L2_zvelc = vnorm_fn(this%jmax, velc_jml_ir)
      
      velc_crit_fn = min( velc_crit_fn, dr / L2_rvelc, r / this%jmax / L2_zvelc )
    end do
    
    deallocate( rvelc_jm_ir , velc_jml_ir )

  end function velc_crit_fn
  
end submodule TimeStepping