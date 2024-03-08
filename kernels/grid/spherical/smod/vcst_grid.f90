submodule (SphericalHarmonics) vcst_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl)                        :: tmp
    real(kind=dbl), pointer               :: gin(:,:,:), gout(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(12, sumNS, grid)
    
    gin(1:6,1:2,1:this%nFourier)  => grid(1:12*this%nFourier)
    gout(1:5,1:2,1:this%nFourier) => grid(1:10*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 2
        tmp = gin(1,i2,i)
        
        gout(1,i2,i) = gin(2,i2,i) * tmp
        gout(2,i2,i) = gin(3,i2,i) * tmp
        gout(3,i2,i) = gin(4,i2,i) * tmp
        gout(4,i2,i) = gin(5,i2,i) * tmp
        gout(5,i2,i) = gin(6,i2,i) * tmp
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(10, grid, sumNS)
    
  end subroutine grid_op_2_vcst_sub
  
  module pure subroutine grid_op_4_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl)                        :: tmp
    real(kind=dbl), pointer               :: gin(:,:,:), gout(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(24, sumNS, grid)
    
    gin(1:6,1:4,1:this%nFourier)  => grid(1:24*this%nFourier)
    gout(1:5,1:4,1:this%nFourier) => grid(1:20*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 4
        tmp = gin(1,i2,i)
        
        gout(1,i2,i) = gin(2,i2,i) * tmp
        gout(2,i2,i) = gin(3,i2,i) * tmp
        gout(3,i2,i) = gin(4,i2,i) * tmp
        gout(4,i2,i) = gin(5,i2,i) * tmp
        gout(5,i2,i) = gin(6,i2,i) * tmp
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(20, grid, sumNS)
    
  end subroutine grid_op_4_vcst_sub
  
  module pure subroutine grid_op_8_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl)                        :: tmp
    real(kind=dbl), pointer               :: gin(:,:,:), gout(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(48, sumNS, grid)
    
    gin(1:6,1:8,1:this%nFourier)  => grid(1:48*this%nFourier)
    gout(1:5,1:8,1:this%nFourier) => grid(1:40*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 8
        tmp = gin(1,i2,i)
        
        gout(1,i2,i) = gin(2,i2,i) * tmp
        gout(2,i2,i) = gin(3,i2,i) * tmp
        gout(3,i2,i) = gin(4,i2,i) * tmp
        gout(4,i2,i) = gin(5,i2,i) * tmp
        gout(5,i2,i) = gin(6,i2,i) * tmp
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(40, grid, sumNS)
    
  end subroutine grid_op_8_vcst_sub
  
  module pure subroutine grid_op_16_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl)                        :: tmp
    real(kind=dbl), pointer               :: gin(:,:,:), gout(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(96, sumNS, grid)
    
    gin(1:6,1:16,1:this%nFourier)  => grid(1:96*this%nFourier)
    gout(1:5,1:16,1:this%nFourier) => grid(1:80*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 16
        tmp = gin(1,i2,i)
        
        gout(1,i2,i) = gin(2,i2,i) * tmp
        gout(2,i2,i) = gin(3,i2,i) * tmp
        gout(3,i2,i) = gin(4,i2,i) * tmp
        gout(4,i2,i) = gin(5,i2,i) * tmp
        gout(5,i2,i) = gin(6,i2,i) * tmp
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(80, grid, sumNS)
    
  end subroutine grid_op_16_vcst_sub
  
end submodule vcst_grid