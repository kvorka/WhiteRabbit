submodule(SphericalHarmonics) vcvv_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(12, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+12*i:12+12*i)
      gout => grid(1+ 2*i: 2+ 2*i)
      
      gout(1) = gin(1) * gin( 4) + gin(2) * gin( 5) + gin(3) * gin( 6)
      gout(2) = gin(7) * gin(10) + gin(8) * gin(11) + gin(9) * gin(12)
    end do
    
    call this%fourtrans%exec_r2c_sub(2, grid, sumNS)
    
  end subroutine grid_op_2_vcvv_sub
  
  module pure subroutine grid_op_4_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(24, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+24*i:24+24*i)
      gout => grid(1+ 4*i: 4+ 4*i)
      
      gout(1) = gin( 1) * gin( 4) + gin( 2) * gin( 5) + gin( 3) * gin( 6)
      gout(2) = gin( 7) * gin(10) + gin( 8) * gin(11) + gin( 9) * gin(12)
      gout(3) = gin(13) * gin(16) + gin(14) * gin(17) + gin(15) * gin(18)
      gout(4) = gin(19) * gin(22) + gin(20) * gin(23) + gin(21) * gin(24)
    end do
    
    call this%fourtrans%exec_r2c_sub(4, grid, sumNS)
    
  end subroutine grid_op_4_vcvv_sub
  
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
  
  module pure subroutine grid_op_16_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(96, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+96*i:96+96*i)
      gout => grid(1+16*i:16+16*i)
      
      gout( 1) = gin( 1) * gin( 4) + gin( 2) * gin( 5) + gin( 3) * gin( 6)
      gout( 2) = gin( 7) * gin(10) + gin( 8) * gin(11) + gin( 9) * gin(12)
      gout( 3) = gin(13) * gin(16) + gin(14) * gin(17) + gin(15) * gin(18)
      gout( 4) = gin(19) * gin(22) + gin(20) * gin(23) + gin(21) * gin(24)
      gout( 5) = gin(25) * gin(28) + gin(26) * gin(29) + gin(27) * gin(30)
      gout( 6) = gin(31) * gin(34) + gin(32) * gin(35) + gin(33) * gin(36)
      gout( 7) = gin(37) * gin(40) + gin(38) * gin(41) + gin(39) * gin(42)
      gout( 8) = gin(43) * gin(46) + gin(44) * gin(47) + gin(45) * gin(48)
      gout( 9) = gin(49) * gin(52) + gin(50) * gin(53) + gin(51) * gin(54)
      gout(10) = gin(55) * gin(58) + gin(56) * gin(59) + gin(57) * gin(60)
      gout(11) = gin(61) * gin(64) + gin(62) * gin(65) + gin(63) * gin(66)
      gout(12) = gin(67) * gin(70) + gin(68) * gin(71) + gin(69) * gin(72)
      gout(13) = gin(73) * gin(76) + gin(74) * gin(77) + gin(75) * gin(78)
      gout(14) = gin(79) * gin(82) + gin(80) * gin(83) + gin(81) * gin(84)
      gout(15) = gin(85) * gin(88) + gin(86) * gin(89) + gin(87) * gin(90)
      gout(16) = gin(91) * gin(94) + gin(92) * gin(95) + gin(93) * gin(96)
    end do
    
    call this%fourtrans%exec_r2c_sub(16, grid, sumNS)
    
  end subroutine grid_op_16_vcvv_sub
  
end submodule vcvv_grid