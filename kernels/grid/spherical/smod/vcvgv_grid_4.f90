submodule (SphericalHarmonics) vcvgv_grid_4
  implicit none; contains
  
  module pure subroutine grid_op_4_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:,:)
    real(kind=dbl), allocatable           :: tmp(:)
    
    call this%fourtrans%exec_c2r_sub(48, sumNS, grid)
    
    gin(1:3,1:4,1:4,1:this%nFourier) => grid(1:48*this%nFourier)
    gout(1:3,1:4,1:this%nFourier)    => grid(1:12*this%nFourier)
    
    allocate( tmp(3) )
    
    do i = 1, this%nFourier
      do i2 = 1, 4
        tmp = gin(1:3,1,i2,i)
        
        gout(1,i2,i) = gin(1,2,i2,i) * tmp(1) + gin(2,2,i2,i) * tmp(2) + gin(3,2,i2,i) * tmp(3)
        gout(2,i2,i) = gin(1,3,i2,i) * tmp(1) + gin(2,3,i2,i) * tmp(2) + gin(3,3,i2,i) * tmp(3)
        gout(3,i2,i) = gin(1,4,i2,i) * tmp(1) + gin(2,4,i2,i) * tmp(2) + gin(3,4,i2,i) * tmp(3)
      end do
    end do
    
    deallocate( tmp )
    
    call this%fourtrans%exec_r2c_sub(12, grid, sumNS)
    
  end subroutine grid_op_4_vcvgv_sub
  
end submodule vcvgv_grid_4