submodule(IceMod) LateralConductivityIce
  implicit none ; contains
  
  module subroutine lambda_ice_jm_sub(this)
    class(T_ice),       intent(inout) :: this
    integer                           :: ir, i1, i2, i3
    real(kind=dbl),     allocatable   :: grid(:,:,:)
    complex(kind=dbl),  allocatable   :: temperature_jm(:), cc_mj(:)
    
    allocate( temperature_jm(this%jms), cc_mj(this%jms2), grid(this%lat_grid%nLegendre,this%lat_grid%nFourier,2) )
      
      !$omp parallel do private(i1, i2, i3, temperature_jm, cc_mj, grid)
      do ir = 1, this%nd
        !Get the temperature field in jm indexing
        call this%temperature_ice_r_jm_sub(ir, temperature_jm)
        
        !Zero the mj indexing holder and reindex the temperature field
        call zero_carray_sub(this%jms2, cc_mj)
        call this%lat_grid%reindexing%scal2scal_jm_to_mj_sub( temperature_jm(1), cc_mj(1), 1, 1 )
        
        !Space to grid :: temperature
        call this%lat_grid%space_to_grid_sub( cc_mj(1), grid )
        
        !Temperature to dimensional units
        grid = (this%Td - this%Tu) * grid + this%Tu
        
        !Compute conductivity on the grid
        do concurrent ( i3 = 1:2, i2 = 1:this%lat_grid%nFourier, i1 = 1:this%lat_grid%nLegendre )
          grid(i1,i2,i3) = name_conductivity_fn( grid(i1,i2,i3) ) / this%lambdaU
        end do
        
        !Grid to space :: conductivity
        call zero_carray_sub(this%jms2, cc_mj)
        call this%lat_grid%grid_to_space_sub( grid, cc_mj(1) )
        
        !Get the viscosity field in jm indexing
        call this%lat_grid%reindexing%scal2scal_mj_to_jm_sub( cc_mj(1), 1, 1, this%mparams%lambda(1,ir), 1, 1)
      end do
      !$omp end parallel do
    
    deallocate( temperature_jm, cc_mj, grid )
    
  end subroutine lambda_ice_jm_sub
  
end submodule LateralConductivityIce