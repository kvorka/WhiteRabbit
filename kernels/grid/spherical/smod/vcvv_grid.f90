submodule(SphericalHarmonics) vcvv_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(12, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 1
        grid(1+(i2+2*i1)) = grid(1+6*(i2+2*i1)) * grid(4+6*(i2+2*i1)) + &
                          & grid(2+6*(i2+2*i1)) * grid(5+6*(i2+2*i1)) + &
                          & grid(3+6*(i2+2*i1)) * grid(6+6*(i2+2*i1))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(2, grid, sumNS)
    
  end subroutine grid_op_2_vcvv_sub
  
  module pure subroutine grid_op_4_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(24, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 3
        grid(1+(i2+4*i1)) = grid(1+6*(i2+4*i1)) * grid(4+6*(i2+4*i1)) + &
                          & grid(2+6*(i2+4*i1)) * grid(5+6*(i2+4*i1)) + &
                          & grid(3+6*(i2+4*i1)) * grid(6+6*(i2+4*i1))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(4, grid, sumNS)
    
  end subroutine grid_op_4_vcvv_sub
  
  module pure subroutine grid_op_8_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(48, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 7
        grid(1+(i2+8*i1)) = grid(1+6*(i2+8*i1)) * grid(4+6*(i2+8*i1)) + &
                          & grid(2+6*(i2+8*i1)) * grid(5+6*(i2+8*i1)) + &
                          & grid(3+6*(i2+8*i1)) * grid(6+6*(i2+8*i1))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(8, grid, sumNS)
    
  end subroutine grid_op_8_vcvv_sub
  
  module pure subroutine grid_op_16_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(96, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 15
        grid(1+(i2+16*i1)) = grid(1+6*(i2+16*i1)) * grid(4+6*(i2+16*i1)) + &
                           & grid(2+6*(i2+16*i1)) * grid(5+6*(i2+16*i1)) + &
                           & grid(3+6*(i2+16*i1)) * grid(6+6*(i2+16*i1))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(16, grid, sumNS)
    
  end subroutine grid_op_16_vcvv_sub
  
end submodule vcvv_grid