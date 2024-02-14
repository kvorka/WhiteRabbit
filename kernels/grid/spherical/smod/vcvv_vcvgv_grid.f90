submodule(SphericalHarmonics) vcvv_vcvgv_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(30, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 1
        grid(1+4*(i2+i1*2)) = grid(1+15*(i2+i1*2)) * grid( 4+15*(i2+i1*2)) + &
                            & grid(2+15*(i2+i1*2)) * grid( 5+15*(i2+i1*2)) + &
                            & grid(3+15*(i2+i1*2)) * grid( 6+15*(i2+i1*2))
        
        grid(2+4*(i2+i1*2)) = grid(4+15*(i2+i1*2)) * grid( 7+15*(i2+i1*2)) + &
                            & grid(5+15*(i2+i1*2)) * grid( 8+15*(i2+i1*2)) + &
                            & grid(6+15*(i2+i1*2)) * grid( 9+15*(i2+i1*2))
        
        grid(3+4*(i2+i1*2)) = grid(4+15*(i2+i1*2)) * grid(10+15*(i2+i1*2)) + &
                            & grid(5+15*(i2+i1*2)) * grid(11+15*(i2+i1*2)) + &
                            & grid(6+15*(i2+i1*2)) * grid(12+15*(i2+i1*2))
        
        grid(4+4*(i2+i1*2)) = grid(4+15*(i2+i1*2)) * grid(13+15*(i2+i1*2)) + &
                            & grid(5+15*(i2+i1*2)) * grid(14+15*(i2+i1*2)) + &
                            & grid(6+15*(i2+i1*2)) * grid(15+15*(i2+i1*2))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(8, this%maxj, grid, sumNS)
    
  end subroutine grid_op_2_vcvv_vcvgv_sub
  
  module pure subroutine grid_op_4_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(60, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 3
        grid(1+4*(i2+i1*4)) = grid(1+15*(i2+i1*4)) * grid( 4+15*(i2+i1*4)) + &
                            & grid(2+15*(i2+i1*4)) * grid( 5+15*(i2+i1*4)) + &
                            & grid(3+15*(i2+i1*4)) * grid( 6+15*(i2+i1*4))
        
        grid(2+4*(i2+i1*4)) = grid(4+15*(i2+i1*4)) * grid( 7+15*(i2+i1*4)) + &
                            & grid(5+15*(i2+i1*4)) * grid( 8+15*(i2+i1*4)) + &
                            & grid(6+15*(i2+i1*4)) * grid( 9+15*(i2+i1*4))
        
        grid(3+4*(i2+i1*4)) = grid(4+15*(i2+i1*4)) * grid(10+15*(i2+i1*4)) + &
                            & grid(5+15*(i2+i1*4)) * grid(11+15*(i2+i1*4)) + &
                            & grid(6+15*(i2+i1*4)) * grid(12+15*(i2+i1*4))
        
        grid(4+4*(i2+i1*4)) = grid(4+15*(i2+i1*4)) * grid(13+15*(i2+i1*4)) + &
                            & grid(5+15*(i2+i1*4)) * grid(14+15*(i2+i1*4)) + &
                            & grid(6+15*(i2+i1*4)) * grid(15+15*(i2+i1*4))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(16, this%maxj, grid, sumNS)
    
  end subroutine grid_op_4_vcvv_vcvgv_sub
  
  module pure subroutine grid_op_8_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(120, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 7
        grid(1+4*(i2+i1*8)) = grid(1+15*(i2+i1*8)) * grid( 4+15*(i2+i1*8)) + &
                            & grid(2+15*(i2+i1*8)) * grid( 5+15*(i2+i1*8)) + &
                            & grid(3+15*(i2+i1*8)) * grid( 6+15*(i2+i1*8))
        
        grid(2+4*(i2+i1*8)) = grid(4+15*(i2+i1*8)) * grid( 7+15*(i2+i1*8)) + &
                            & grid(5+15*(i2+i1*8)) * grid( 8+15*(i2+i1*8)) + &
                            & grid(6+15*(i2+i1*8)) * grid( 9+15*(i2+i1*8))
        
        grid(3+4*(i2+i1*8)) = grid(4+15*(i2+i1*8)) * grid(10+15*(i2+i1*8)) + &
                            & grid(5+15*(i2+i1*8)) * grid(11+15*(i2+i1*8)) + &
                            & grid(6+15*(i2+i1*8)) * grid(12+15*(i2+i1*8))
        
        grid(4+4*(i2+i1*8)) = grid(4+15*(i2+i1*8)) * grid(13+15*(i2+i1*8)) + &
                            & grid(5+15*(i2+i1*8)) * grid(14+15*(i2+i1*8)) + &
                            & grid(6+15*(i2+i1*8)) * grid(15+15*(i2+i1*8))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(32, this%maxj, grid, sumNS)
    
  end subroutine grid_op_8_vcvv_vcvgv_sub
  
  module pure subroutine grid_op_16_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i1, i2
    
    call this%fourtrans%exec_c2r_sub(240, this%maxj+1, sumNS, grid)
    
    do i1 = 0, this%nFourier-1
      do i2 = 0, 15
        grid(1+4*(i2+i1*16)) = grid(1+15*(i2+i1*16)) * grid( 4+15*(i2+i1*16)) + &
                             & grid(2+15*(i2+i1*16)) * grid( 5+15*(i2+i1*16)) + &
                             & grid(3+15*(i2+i1*16)) * grid( 6+15*(i2+i1*16))
        
        grid(2+4*(i2+i1*16)) = grid(4+15*(i2+i1*16)) * grid( 7+15*(i2+i1*16)) + &
                             & grid(5+15*(i2+i1*16)) * grid( 8+15*(i2+i1*16)) + &
                             & grid(6+15*(i2+i1*16)) * grid( 9+15*(i2+i1*16))
        
        grid(3+4*(i2+i1*16)) = grid(4+15*(i2+i1*16)) * grid(10+15*(i2+i1*16)) + &
                             & grid(5+15*(i2+i1*16)) * grid(11+15*(i2+i1*16)) + &
                             & grid(6+15*(i2+i1*16)) * grid(12+15*(i2+i1*16))
        
        grid(4+4*(i2+i1*16)) = grid(4+15*(i2+i1*16)) * grid(13+15*(i2+i1*16)) + &
                             & grid(5+15*(i2+i1*16)) * grid(14+15*(i2+i1*16)) + &
                             & grid(6+15*(i2+i1*16)) * grid(15+15*(i2+i1*16))
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(64, this%maxj, grid, sumNS)
    
  end subroutine grid_op_16_vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv_grid