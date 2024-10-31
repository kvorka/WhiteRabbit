submodule (physicalobject) meassures
  implicit none ; contains
  
  module procedure nuss_fn
    
    nuss_fn = c2r_fn( -this%q_r_fn(this%nd,1,1) ) / ( this%r_ud * s4pi )
    
  end procedure nuss_fn
  
  module procedure reynolds_fn
    integer                        :: ir
    real(kind=dbl),    allocatable :: field_vals(:)
    complex(kind=dbl), allocatable :: velocity(:)
    
    allocate( field_vals(this%nd+1), velocity(this%jmv) )
      
      if ( ( present(choice) ) .and. ( choice == 'convective' ) ) then
        
        !$omp parallel do private(velocity)
        do ir = 1, this%nd+1
          call this%sol%conv_velocity_jml_sub( ir, velocity )
          field_vals(ir) = vectnorm2_fn( this%jmax, velocity )
        end do
        !$omp end parallel do
        
      else
        
        !$omp parallel do private(velocity)
        do ir = 1, this%nd+1
          call this%v_rr_ijml_sub( ir, velocity )
          field_vals(ir) = vectnorm2_fn( this%jmax, velocity )
        end do
        !$omp end parallel do
        
      end if
      
      reynolds_fn = sqrt( this%rad_grid%intV_fn( field_vals ) / this%rad_grid%volume )
      
    deallocate( field_vals, velocity )
    
  end procedure reynolds_fn

end submodule meassures