submodule(IceMod) LateralViscosityIce
  implicit none ; contains
  
  module subroutine visc_ice_jm_sub(this)
    class(T_ice),       intent(inout) :: this
    integer                           :: ir, i1, i2, i3
    real(kind=dbl)                    :: avrgstress
    real(kind=dbl),     allocatable   :: grid(:,:,:)
    complex(kind=dbl),  allocatable   :: temperature_jm(:), cc_mj(:)
    
    allocate( temperature_jm(this%jms), cc_mj(this%jms2), grid(this%lat_grid%nLegendre,this%lat_grid%nFourier,2) )
      
      !$omp parallel do private(i1, i2, i3, avrgstress, temperature_jm, cc_mj, grid)
      do ir = 1, this%nd
        !Usefull for viscosity computation
        avrgstress = this%devstress_ice_r_fn(ir)
        
        !Get the temperature field in jm indexing
        call this%temp_ice_r_jm_sub(ir, temperature_jm)
        
        !Zero the mj indexing holder and reindex the temperature field
        call zero_carray_sub(this%jms2, cc_mj)
        call this%lat_grid%reindexing%scal2scal_jm_to_mj_sub( temperature_jm(1), cc_mj(1), 1, 1 )
        
        !Space to grid :: temperature
        call this%lat_grid%space_to_grid_sub( cc_mj(1), grid )
        
        !Temperature to dimensional units
        grid = (this%Td - this%Tu) * grid + this%Tu
        
        !Compute viscosity on the grid
        if ( this%andrade ) then
          do concurrent ( i3 = 1:2, i2 = 1:this%lat_grid%nFourier, i1 = 1:this%lat_grid%nLegendre )
            grid(i1,i2,i3) = min( goldsby_visc_fn( this%diam, grid(i1,i2,i3), avrgstress ), this%cutoff )
            grid(i1,i2,i3) = andrade_visc_fn( this%mu, this%omega, grid(i1,i2,i3) ) / this%viscU
          end do
        else
          do concurrent ( i3 = 1:2, i2 = 1:this%lat_grid%nFourier, i1 = 1:this%lat_grid%nLegendre )
            grid(i1,i2,i3) = min( goldsby_visc_fn( this%diam, grid(i1,i2,i3), avrgstress ), this%cutoff )
            grid(i1,i2,i3) = grid(i1,i2,i3) / this%viscU
          end do
        end if
        
        !Grid to space :: viscosity
        call zero_carray_sub(this%jms2, cc_mj)
        call this%lat_grid%grid_to_space_sub( grid, cc_mj(1) )
        
        !Get the viscosity field in jm indexing
        call this%lat_grid%reindexing%scal2scal_mj_to_jm_sub( cc_mj(1), 1, 1, this%mparams%visc(1,ir), 1, 1)
      end do
      !$omp end parallel do
    
    deallocate( temperature_jm, cc_mj, grid )
    
  end subroutine visc_ice_jm_sub
  
end submodule LateralViscosityIce