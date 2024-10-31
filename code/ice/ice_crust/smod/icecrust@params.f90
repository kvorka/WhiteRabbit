submodule (icecrust) params
  implicit none; contains
  
  module procedure lambda_iceCrust_jm_sub
    integer                        :: ir, ijm, i1, i2, i3
    real(kind=dbl),    allocatable :: grid(:,:,:)
    complex(kind=dbl), allocatable :: temp_jm(:), cc_mj(:)
    
    if ( this%mparams%initlambda ) then

      if ( this%mparams%lambda_radial) then
        !$omp parallel do
        do ir = 1, this%nd
          if ( this%rad_grid%r(ir) < this%ru - this%hC ) then
            this%mparams%lambda(1,ir) = r2c_fn( s4pi * this%lambdaU / name_conductivity_fn( this%avrg_temperature_ice_ir_fn(ir) ) )
          else
            this%mparams%lambda(1,ir) = r2c_fn( s4pi * this%lambdaU / this%lambdaC )
          end if
        end do
        !$omp end parallel do
        
      else
        call this%lat_grid%reindexing%allocate_scalars_sub(1, temp_jm)
        call this%lat_grid%reindexing%allocate_scalars_sub(1, cc_mj)
        call this%lat_grid%allocate_grid_sub( grid )
        
        !$omp parallel do private(ijm, i1, i2, i3, temp_jm, cc_mj, grid)
        do ir = 1, this%nd
          !Get the non-dimensional temperature field in jm indexing
          call this%temp_r_ijm_sub(ir, temp_jm)
          
          !Zero the mj indexing holder and reindex the temperature field
          call zero_carray_sub( this%jms2, cc_mj(1) )
          call this%lat_grid%reindexing%scal2scal_jm_to_mj_sub( temp_jm(1), cc_mj(1), 1, 1 )
          
          !Space to grid :: non-dimensional temperature
          call this%lat_grid%space_to_grid_sub( cc_mj(1), grid )
          
          !Compute non-dimensional 1/conductivity on the grid
          do concurrent ( i3 = 1:2, i2 = 1:size(grid,dim=2), i1 = 1:size(grid,dim=1) )
            grid(i1,i2,i3) = ( this%Td - this%Tu ) * grid(i1,i2,i3) + this%Tu
            grid(i1,i2,i3) = this%lambdaU / name_conductivity_fn( grid(i1,i2,i3) )
          end do
          
          !Grid to space :: non-dimensional conductivity
          call zero_carray_sub(this%jms2, cc_mj)
          call this%lat_grid%grid_to_space_sub( grid, cc_mj(1) )
          
          !Get the non-dimensional conductivity field in jm indexing
          call this%lat_grid%reindexing%scal2scal_mj_to_jm_sub( cc_mj(1), 1, 1, this%mparams%lambda(1,ir), 1, 1)
        end do
        !$omp end parallel do
        
        deallocate( temp_jm, cc_mj, grid )
        
      end if
    end if
    
  end procedure lambda_iceCrust_jm_sub
  
  module procedure cp_iceCrust_jm_sub
    integer                        :: ir, ijm, i1, i2, i3
    real(kind=dbl),    allocatable :: grid(:,:,:)
    complex(kind=dbl), allocatable :: temp_jm(:), cc_mj(:)
    
    if ( this%mparams%initcp ) then
      
      if ( this%mparams%cp_radial) then
        !$omp parallel do
        do ir = 1, this%nd+1
          this%mparams%cp(1,ir) = r2c_fn( s4pi * this%cU / name_cp_fn( this%avrg_temperature_ice_irr_fn(ir) ) )
        end do
        !$omp end parallel do
        
      else
        call this%lat_grid%reindexing%allocate_scalars_sub(1, temp_jm)
        call this%lat_grid%reindexing%allocate_scalars_sub(1, cc_mj)
        call this%lat_grid%allocate_grid_sub( grid )
        
        !$omp parallel do private(ijm, i1, i2, i3, temp_jm, cc_mj, grid)
        do ir = 1, this%nd+1
          !Get the non-dimensional temperature field in jm indexing
          call this%temp_rr_ijm_sub(ir, temp_jm)
          
          !Zero the mj indexing holder and reindex the temperature field
          call zero_carray_sub( this%jms2, cc_mj(1) )
          call this%lat_grid%reindexing%scal2scal_jm_to_mj_sub( temp_jm(1), cc_mj(1), 1, 1 )
          
          !Space to grid :: non-dimensional temperature
          call this%lat_grid%space_to_grid_sub( cc_mj(1), grid )
          
          !Compute non-dimensional 1/capacity on the grid
          do concurrent ( i3 = 1:2, i2 = 1:size(grid,dim=2), i1 = 1:size(grid,dim=1) )
            grid(i1,i2,i3) = ( this%Td - this%Tu ) * grid(i1,i2,i3) + this%Tu
            grid(i1,i2,i3) = this%cU / name_cp_fn( grid(i1,i2,i3) )
          end do
          
          !Grid to space :: non-dimensional 1/capacity
          call zero_carray_sub(this%jms2, cc_mj)
          call this%lat_grid%grid_to_space_sub( grid, cc_mj(1) )
          
          !Get the non-dimensional 1/capacity field in jm indexing
          call this%lat_grid%reindexing%scal2scal_mj_to_jm_sub( cc_mj(1), 1, 1, this%mparams%cp(1,ir), 1, 1)
        end do
        !$omp end parallel do
        
        deallocate( temp_jm, cc_mj, grid )
        
      end if
    end if
    
  end procedure cp_iceCrust_jm_sub
  
  module procedure alpha_iceCrust_jm_sub
    integer                        :: ir, ijm, i1, i2, i3
    real(kind=dbl),    allocatable :: grid(:,:,:)
    complex(kind=dbl), allocatable :: temp_jm(:), cc_mj(:)
    
    if ( this%mparams%initalpha ) then

      if ( this%mparams%alpha_radial) then
        !$omp parallel do
        do ir = 1, this%nd+1
          this%mparams%alpha(1,ir) = r2c_fn( s4pi * name_expansivity_fn( this%avrg_temperature_ice_irr_fn(ir) ) / this%alphaU )
        end do
        !$omp end parallel do
        
      else
        call this%lat_grid%reindexing%allocate_scalars_sub(1, temp_jm)
        call this%lat_grid%reindexing%allocate_scalars_sub(1, cc_mj)
        call this%lat_grid%allocate_grid_sub( grid )
        
        !$omp parallel do private(ijm, i1, i2, i3, temp_jm, cc_mj, grid)
        do ir = 1, this%nd+1
          !Get the non-dimensional temperature field in jm indexing
          call this%temp_rr_ijm_sub(ir, temp_jm)
          
          !Zero the mj indexing holder and reindex the temperature field
          call zero_carray_sub( this%jms2, cc_mj(1) )
          call this%lat_grid%reindexing%scal2scal_jm_to_mj_sub( temp_jm(1), cc_mj(1), 1, 1 )
          
          !Space to grid :: non-dimensional temperature
          call this%lat_grid%space_to_grid_sub( cc_mj(1), grid )
          
          !Compute non-dimensional expansivity on the grid
          do concurrent ( i3 = 1:2, i2 = 1:size(grid,dim=2), i1 = 1:size(grid,dim=1) )
            grid(i1,i2,i3) = ( this%Td - this%Tu ) * grid(i1,i2,i3) + this%Tu
            grid(i1,i2,i3) = name_expansivity_fn( grid(i1,i2,i3) ) / this%alphaU
          end do
          
          !Grid to space :: non-dimensional conductivity
          call zero_carray_sub(this%jms2, cc_mj)
          call this%lat_grid%grid_to_space_sub( grid, cc_mj(1) )
          
          !Get the non-dimensional conductivity field in jm indexing
          call this%lat_grid%reindexing%scal2scal_mj_to_jm_sub( cc_mj(1), 1, 1, this%mparams%cp(1,ir), 1, 1)
        end do
        !$omp end parallel do
        
        deallocate( temp_jm, cc_mj, grid )
        
      end if
    end if
    
  end procedure alpha_iceCrust_jm_sub
  
  module procedure visc_iceCrust_jm_sub
    integer                        :: ir, i1, i2, i3
    real(kind=dbl)                 :: avrgstress
    real(kind=dbl),    allocatable :: grid(:,:,:)
    complex(kind=dbl), allocatable :: temp_jm(:), cc_mj(:)
    
    if ( this%mparams%initvisc ) then
      if ( this%mparams%visc_radial) then
        !$omp parallel do
        do ir = 1, this%nd
          this%mparams%visc(1,ir) = r2c_fn( s4pi * this%viscU / min( goldsby_visc_fn( this%diam,                           &
                                                                                    & this%avrg_temperature_ice_ir_fn(ir), &
                                                                                    & this%avrg_stress_ice_ir_fn(ir)       ), &
                                                                   & this%cutoff                                              ) )
        end do
        !$omp end parallel do
        
      else
        call this%lat_grid%reindexing%allocate_scalars_sub(1, temp_jm)
        call this%lat_grid%reindexing%allocate_scalars_sub(1, cc_mj)
        call this%lat_grid%allocate_grid_sub( grid )
        
        !$omp parallel do private(i1, i2, i3, avrgstress, temp_jm, cc_mj, grid)
        do ir = 1, this%nd
          !Usefull for viscosity computation
          avrgstress = this%avrg_stress_ice_ir_fn(ir)
          
          !Get the temperature field in jm indexing
          call this%temp_r_ijm_sub(ir, temp_jm)
          
          !Zero the mj indexing holder and reindex the temperature field
          call zero_carray_sub(this%jms2, cc_mj)
          call this%lat_grid%reindexing%scal2scal_jm_to_mj_sub( temp_jm(1), cc_mj(1), 1, 1 )
          
          !Space to grid :: temperature
          call this%lat_grid%space_to_grid_sub( cc_mj(1), grid )
          
          !Compute 1/viscosity on the grid
          do concurrent ( i3 = 1:2, i2 = 1:size(grid,dim=2), i1 = 1:size(grid,dim=1) )
            grid(i1,i2,i3) = ( this%Td - this%Tu ) * grid(i1,i2,i3) + this%Tu
            grid(i1,i2,i3) = this%viscU / min( goldsby_visc_fn( this%diam, grid(i1,i2,i3), avrgstress ), this%cutoff )
          end do
          
          !Grid to space :: viscosity
          call zero_carray_sub(this%jms2, cc_mj)
          call this%lat_grid%grid_to_space_sub( grid, cc_mj(1) )
          
          !Get the viscosity field in jm indexing
          call this%lat_grid%reindexing%scal2scal_mj_to_jm_sub( cc_mj(1), 1, 1, this%mparams%visc(1,ir), 1, 1)
        end do
        !$omp end parallel do
        
        deallocate( temp_jm, cc_mj, grid )

      end if
    end if
    
  end procedure visc_iceCrust_jm_sub
  
  module procedure surfTemp_iceCrust_jm_sub
    integer                        :: i1, i2
    real(kind=dbl)                 :: theta
    real(kind=dbl),    allocatable :: grid(:,:,:)
    complex(kind=dbl), allocatable :: cc_mj(:)
    
    call this%lat_grid%reindexing%allocate_scalars_sub(1, cc_mj)
    call this%lat_grid%allocate_grid_sub( grid )
    
    !to grid :: non-dimensional surface temperature
    do concurrent ( i2 = 1:size(grid,dim=2), i1 = 1:size(grid,dim=1) )
      grid(i1,i2,1) = ( name_surfaceTemp_fn(   acos(this%lat_grid%lgp%rootsweights(1,i1))) - this%Tu ) / ( this%Td - this%Tu )
      grid(i1,i2,2) = ( name_surfaceTemp_fn(pi-acos(this%lat_grid%lgp%rootsweights(1,i1))) - this%Tu ) / ( this%Td - this%Tu )
    end do
    
    !Grid to space :: non-dim surface temperature
    call zero_carray_sub(this%jms2, cc_mj)
    call this%lat_grid%grid_to_space_sub( grid, cc_mj(1) )
      
    !Get the surface temperature field in jm indexing
    call this%lat_grid%reindexing%scal2scal_mj_to_jm_sub( cc_mj(1), 1, 1, this%bnd%temp_up(1), 1, 1)
    
    deallocate( cc_mj, grid )
    
  end procedure surfTemp_iceCrust_jm_sub
  
end submodule params