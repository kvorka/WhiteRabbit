submodule(SphericalHarmonics) vcvgv_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(24, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 1
        grid(1+3*(i2+2*i1)) = grid(1+12*(i2+2*i1)) * grid(10+12*(i2+2*i1)) + &
                            & grid(2+12*(i2+2*i1)) * grid(11+12*(i2+2*i1)) + &
                            & grid(3+12*(i2+2*i1)) * grid(12+12*(i2+2*i1))
        
        grid(2+3*(i2+2*i1)) = grid(4+12*(i2+2*i1)) * grid(10+12*(i2+2*i1)) + &
                            & grid(5+12*(i2+2*i1)) * grid(11+12*(i2+2*i1)) + &
                            & grid(6+12*(i2+2*i1)) * grid(12+12*(i2+2*i1))
        
        grid(3+3*(i2+2*i1)) = grid(7+12*(i2+2*i1)) * grid(10+12*(i2+2*i1)) + &
                            & grid(8+12*(i2+2*i1)) * grid(11+12*(i2+2*i1)) + &
                            & grid(9+12*(i2+2*i1)) * grid(12+12*(i2+2*i1))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(6, this%maxj+1, grid, sumNS)
    
  end subroutine grid_op_2_vcvgv_sub
  
  module pure subroutine grid_op_4_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(48, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 3
        grid(1+3*(i2+4*i1)) = grid(1+12*(i2+4*i1)) * grid(10+12*(i2+4*i1)) + &
                            & grid(2+12*(i2+4*i1)) * grid(11+12*(i2+4*i1)) + &
                            & grid(3+12*(i2+4*i1)) * grid(12+12*(i2+4*i1))
        
        grid(2+3*(i2+4*i1)) = grid(4+12*(i2+4*i1)) * grid(10+12*(i2+4*i1)) + &
                            & grid(5+12*(i2+4*i1)) * grid(11+12*(i2+4*i1)) + &
                            & grid(6+12*(i2+4*i1)) * grid(12+12*(i2+4*i1))
        
        grid(3+3*(i2+4*i1)) = grid(7+12*(i2+4*i1)) * grid(10+12*(i2+4*i1)) + &
                            & grid(8+12*(i2+4*i1)) * grid(11+12*(i2+4*i1)) + &
                            & grid(9+12*(i2+4*i1)) * grid(12+12*(i2+4*i1))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(12, this%maxj+1, grid, sumNS)
    
  end subroutine grid_op_4_vcvgv_sub
  
  module pure subroutine grid_op_8_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(96, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 7
        grid(1+3*(i2+8*i1)) = grid(1+12*(i2+8*i1)) * grid(10+12*(i2+8*i1)) + &
                            & grid(2+12*(i2+8*i1)) * grid(11+12*(i2+8*i1)) + &
                            & grid(3+12*(i2+8*i1)) * grid(12+12*(i2+8*i1))
        
        grid(2+3*(i2+8*i1)) = grid(4+12*(i2+8*i1)) * grid(10+12*(i2+8*i1)) + &
                            & grid(5+12*(i2+8*i1)) * grid(11+12*(i2+8*i1)) + &
                            & grid(6+12*(i2+8*i1)) * grid(12+12*(i2+8*i1))
        
        grid(3+3*(i2+8*i1)) = grid(7+12*(i2+8*i1)) * grid(10+12*(i2+8*i1)) + &
                            & grid(8+12*(i2+8*i1)) * grid(11+12*(i2+8*i1)) + &
                            & grid(9+12*(i2+8*i1)) * grid(12+12*(i2+8*i1))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(24, this%maxj+1, grid, sumNS)
    
  end subroutine grid_op_8_vcvgv_sub
  
  module pure subroutine grid_op_16_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(192, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 15
        grid(1+3*(i2+16*i1)) = grid(1+12*(i2+16*i1)) * grid(10+12*(i2+16*i1)) + &
                             & grid(2+12*(i2+16*i1)) * grid(11+12*(i2+16*i1)) + &
                             & grid(3+12*(i2+16*i1)) * grid(12+12*(i2+16*i1))
        
        grid(2+3*(i2+16*i1)) = grid(4+12*(i2+16*i1)) * grid(10+12*(i2+16*i1)) + &
                             & grid(5+12*(i2+16*i1)) * grid(11+12*(i2+16*i1)) + &
                             & grid(6+12*(i2+16*i1)) * grid(12+12*(i2+16*i1))
        
        grid(3+3*(i2+16*i1)) = grid(7+12*(i2+16*i1)) * grid(10+12*(i2+16*i1)) + &
                             & grid(8+12*(i2+16*i1)) * grid(11+12*(i2+16*i1)) + &
                             & grid(9+12*(i2+16*i1)) * grid(12+12*(i2+16*i1))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(48, this%maxj+1, grid, sumNS)
    
  end subroutine grid_op_16_vcvgv_sub
  
end submodule vcvgv_grid