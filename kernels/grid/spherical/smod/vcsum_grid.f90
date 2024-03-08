submodule(SphericalHarmonics) vcsum_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(4, sumNS, grid)
    
    gin(1:2,1:2,1:this%nFourier) => grid(1:4*this%nFourier)
    gout(1:2,1:this%nFourier)    => grid(1:2*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 2
        gout(i2,i) = gin(1,i2,i) * gin(2,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(2, grid, sumNS)
    
  end subroutine grid_op_2_vcsum_sub
  
  module pure subroutine grid_op_4_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(8, sumNS, grid)
    
    gin(1:2,1:4,1:this%nFourier) => grid(1:8*this%nFourier)
    gout(1:4,1:this%nFourier)    => grid(1:4*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 4
        gout(i2,i) = gin(1,i2,i) * gin(2,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(4, grid, sumNS)
    
  end subroutine grid_op_4_vcsum_sub
  
  module pure subroutine grid_op_8_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(16, sumNS, grid)
    
    gin(1:2,1:8,1:this%nFourier) => grid(1:16*this%nFourier)
    gout(1:8,1:this%nFourier)    => grid(1:8*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 8
        gout(i2,i) = gin(1,i2,i) * gin(2,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(8, grid, sumNS)
    
  end subroutine grid_op_8_vcsum_sub
  
  module pure subroutine grid_op_16_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(32, sumNS, grid)
    
    gin(1:2,1:16,1:this%nFourier) => grid(1:32*this%nFourier)
    gout(1:16,1:this%nFourier)    => grid(1:16*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 16
        gout(i2,i) = gin(1,i2,i) * gin(2,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(16, grid, sumNS)
    
  end subroutine grid_op_16_vcsum_sub
  
end submodule vcsum_grid