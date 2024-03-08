submodule (SphericalHarmonics) vcvv_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(12, sumNS, grid)
    
    gin(1:6,1:2,1:this%nFourier) => grid(1:12*this%nFourier)
    gout(1:2,1:this%nFourier)   => grid(1:2*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 2
        gout(i2,i) = gin(1,i2,i) * gin(4,i2,i) + gin(2,i2,i) * gin(5,i2,i) + gin(3,i2,i) * gin(6,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(2, grid, sumNS)
    
  end subroutine grid_op_2_vcvv_sub
  
  module pure subroutine grid_op_4_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(24, sumNS, grid)
    
    gin(1:6,1:4,1:this%nFourier) => grid(1:24*this%nFourier)
    gout(1:4,1:this%nFourier)   => grid(1:4*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 4
        gout(i2,i) = gin(1,i2,i) * gin(4,i2,i) + gin(2,i2,i) * gin(5,i2,i) + gin(3,i2,i) * gin(6,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(4, grid, sumNS)
    
  end subroutine grid_op_4_vcvv_sub
  
  module pure subroutine grid_op_8_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(48, sumNS, grid)
    
    gin(1:6,1:8,1:this%nFourier) => grid(1:48*this%nFourier)
    gout(1:8,1:this%nFourier)   => grid(1:8*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 8
        gout(i2,i) = gin(1,i2,i) * gin(4,i2,i) + gin(2,i2,i) * gin(5,i2,i) + gin(3,i2,i) * gin(6,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(8, grid, sumNS)
    
  end subroutine grid_op_8_vcvv_sub
  
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
  
end submodule vcvv_grid