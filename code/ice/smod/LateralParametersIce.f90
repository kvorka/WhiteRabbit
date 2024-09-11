submodule(IceMod) LateralParametersIce
  implicit none ; contains
  
  pure subroutine temp_to_diffvisc_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:)
    
    gin(1:step,1:nfour)  => grid(1:step*nfour)
    gout(1:step,1:nfour) => grid(1:step*nfour)
    
    do i1 = 1, nfour
      do i2 = 1, step
        gout(i2,i1) = goldsby_diffvisc_fn(diam_ice, gin(i2,i1))
      end do
    end do
    
  end subroutine temp_to_diffvisc_sub
  
  module subroutine visc_ice_jm_sub(this)
    class(T_ice),      intent(inout) :: this
    integer                          :: ir
    complex(kind=dbl), allocatable   :: temperature(:), viscosity(:)
    
    allocate( temperature(this%jms2), viscosity(this%jms2) )
      
      !$omp parallel do private(temperature, viscosity)
      do ir = 1, this%nd
        call this%temperature_ice_r_jm_sub( ir, temperature )
        
        call this%lat_grid%transform_sub(1, 1, temperature, viscosity, temp_to_diffvisc_sub )
        
        this%sol%visc(:,ir) = viscosity(1:this%jms)
      end do
      !$omp end parallel do
    
    deallocate( temperature, viscosity )
    
  end subroutine visc_ice_jm_sub
  
end submodule LateralParametersIce