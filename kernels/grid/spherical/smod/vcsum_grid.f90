submodule(SphericalHarmonics) vcsum_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(4, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 1
        grid(1+(i2+i1*2)) = grid(1+2*(i2+i1*2)) * grid(2+2*(i2+i1*2))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(2, this%maxj, grid, sumNS)
    
  end subroutine grid_op_2_vcsum_sub
  
  module pure subroutine grid_op_4_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(8, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 3
        grid(1+(i2+i1*4)) = grid(1+2*(i2+i1*4)) * grid(2+2*(i2+i1*4))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(4, this%maxj, grid, sumNS)
    
  end subroutine grid_op_4_vcsum_sub
  
  module pure subroutine grid_op_8_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(16, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 7
        grid(1+(i2+i1*8)) = grid(1+2*(i2+i1*8)) * grid(2+2*(i2+i1*8))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(8, this%maxj, grid, sumNS)
    
  end subroutine grid_op_8_vcsum_sub
  
  module pure subroutine grid_op_16_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(32, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 15
        grid(1+(i2+i1*16)) = grid(1+2*(i2+i1*16)) * grid(2+2*(i2+i1*16))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(16, this%maxj, grid, sumNS)
    
  end subroutine grid_op_16_vcsum_sub
  
end submodule vcsum_grid