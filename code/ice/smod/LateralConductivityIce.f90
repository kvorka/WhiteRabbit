submodule(IceMod) LateralConductivityIce
  implicit none ; contains
  
  module subroutine lambda_ice_jm_sub(this)
    class(T_ice),       intent(inout) :: this
    integer                           :: ir, ijm, i1, i2, i3
    real(kind=dbl),     allocatable   :: grid(:,:,:)
    complex(kind=dbl),  allocatable   :: temp_jm(:), cc_mj(:)
    
    allocate( temp_jm(this%jms), cc_mj(this%jms2), grid(this%lat_grid%nLegendre,this%lat_grid%nFourier,2) )
      
      !$omp parallel do private(ijm, i1, i2, i3, temp_jm, cc_mj, grid)
      do ir = 1, this%nd
        !Get the non-dimensional temperature field in jm indexing
        do concurrent ( ijm = 1:this%jms )
          temp_jm(ijm) = this%temp_r_fn(ir,ijm)
        end do
        
        !Zero the mj indexing holder and reindex the temperature field
        call zero_carray_sub( this%jms2, cc_mj(1) )
        call this%lat_grid%reindexing%scal2scal_jm_to_mj_sub( temp_jm(1), cc_mj(1), 1, 1 )
        
        !Space to grid :: non-dimensional temperature
        call this%lat_grid%space_to_grid_sub( cc_mj(1), grid )
        
        !Compute non-dimensional conductivity on the grid
        do concurrent ( i3 = 1:2, i2 = 1:this%lat_grid%nFourier, i1 = 1:this%lat_grid%nLegendre )
          grid(i1,i2,i3) = ( this%Td - this%Tu ) * grid(i1,i2,i3) + this%Tu
          grid(i1,i2,i3) = name_conductivity_fn( grid(i1,i2,i3) ) / this%lambdaU
        end do
        
        !Grid to space :: non-dimensional conductivity
        call zero_carray_sub(this%jms2, cc_mj)
        call this%lat_grid%grid_to_space_sub( grid, cc_mj(1) )
        
        !Get the non-dimensional conductivity field in jm indexing
        call this%lat_grid%reindexing%scal2scal_mj_to_jm_sub( cc_mj(1), 1, 1, this%mparams%lambda(1,ir), 1, 1)
      end do
      !$omp end parallel do
    
    deallocate( temp_jm, cc_mj, grid )
    
  end subroutine lambda_ice_jm_sub
  
end submodule LateralConductivityIce