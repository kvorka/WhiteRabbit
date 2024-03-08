submodule (SphericalHarmonics) vcvv_grid_16
  implicit none; contains
  
  module pure subroutine grid_op_16_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(96, sumNS, grid)
    
    gin(1:6,1:16,1:this%nFourier) => grid(1:96*this%nFourier)
    gout(1:16,1:this%nFourier)   => grid(1:16*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 16
        gout(i2,i) = gin(1,i2,i) * gin(4,i2,i) + gin(2,i2,i) * gin(5,i2,i) + gin(3,i2,i) * gin(6,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(16, grid, sumNS)
    
  end subroutine grid_op_16_vcvv_sub
  
end submodule vcvv_grid_16