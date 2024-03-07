submodule (SphericalHarmonics) vcvv_grid_8
  implicit none; contains
  
  module pure subroutine grid_op_8_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(48, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+48*i:48+48*i)
      gout => grid(1+ 8*i: 8+ 8*i)
      
      gout(1) = gin( 1) * gin( 4) + gin( 2) * gin( 5) + gin( 3) * gin( 6)
      gout(2) = gin( 7) * gin(10) + gin( 8) * gin(11) + gin( 9) * gin(12)
      gout(3) = gin(13) * gin(16) + gin(14) * gin(17) + gin(15) * gin(18)
      gout(4) = gin(19) * gin(22) + gin(20) * gin(23) + gin(21) * gin(24)
      gout(5) = gin(25) * gin(28) + gin(26) * gin(29) + gin(27) * gin(30)
      gout(6) = gin(31) * gin(34) + gin(32) * gin(35) + gin(33) * gin(36)
      gout(7) = gin(37) * gin(40) + gin(38) * gin(41) + gin(39) * gin(42)
      gout(8) = gin(43) * gin(46) + gin(44) * gin(47) + gin(45) * gin(48)
    end do
    
    call this%fourtrans%exec_r2c_sub(8, grid, sumNS)
    
  end subroutine grid_op_8_vcvv_sub
  
end submodule vcvv_grid_8