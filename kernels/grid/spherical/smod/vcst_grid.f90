submodule (SphericalHarmonics) vcst_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(12, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+12*i:12+12*i)
      gout => grid(1+10*i:10+10*i)
      
      gout( 1) = gin( 1) * gin( 6)
      gout( 2) = gin( 2) * gin( 6)
      gout( 3) = gin( 3) * gin( 6)
      gout( 4) = gin( 4) * gin( 6)
      gout( 5) = gin( 5) * gin( 6)
      gout( 6) = gin( 7) * gin(12)
      gout( 7) = gin( 8) * gin(12)
      gout( 8) = gin( 9) * gin(12)
      gout( 9) = gin(10) * gin(12)
      gout(10) = gin(11) * gin(12)
    end do
    
    call this%fourtrans%exec_r2c_sub(10, grid, sumNS)
    
  end subroutine grid_op_2_vcst_sub
  
  module pure subroutine grid_op_4_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(24, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+24*i:24+24*i)
      gout => grid(1+20*i:20+20*i)
      
      gout( 1) = gin( 1) * gin( 6)
      gout( 2) = gin( 2) * gin( 6)
      gout( 3) = gin( 3) * gin( 6)
      gout( 4) = gin( 4) * gin( 6)
      gout( 5) = gin( 5) * gin( 6)
      gout( 6) = gin( 7) * gin(12)
      gout( 7) = gin( 8) * gin(12)
      gout( 8) = gin( 9) * gin(12)
      gout( 9) = gin(10) * gin(12)
      gout(10) = gin(11) * gin(12)
      gout(11) = gin(13) * gin(18)
      gout(12) = gin(14) * gin(18)
      gout(13) = gin(15) * gin(18)
      gout(14) = gin(16) * gin(18)
      gout(15) = gin(17) * gin(18)
      gout(16) = gin(19) * gin(24)
      gout(17) = gin(20) * gin(24)
      gout(18) = gin(21) * gin(24)
      gout(19) = gin(22) * gin(24)
      gout(20) = gin(23) * gin(24)
    end do
    
    call this%fourtrans%exec_r2c_sub(20, grid, sumNS)
    
  end subroutine grid_op_4_vcst_sub
  
  module pure subroutine grid_op_8_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(48, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+48*i:48+48*i)
      gout => grid(1+40*i:40+40*i)
      
      gout( 1) = gin( 1) * gin( 6)
      gout( 2) = gin( 2) * gin( 6)
      gout( 3) = gin( 3) * gin( 6)
      gout( 4) = gin( 4) * gin( 6)
      gout( 5) = gin( 5) * gin( 6)
      gout( 6) = gin( 7) * gin(12)
      gout( 7) = gin( 8) * gin(12)
      gout( 8) = gin( 9) * gin(12)
      gout( 9) = gin(10) * gin(12)
      gout(10) = gin(11) * gin(12)
      gout(11) = gin(13) * gin(18)
      gout(12) = gin(14) * gin(18)
      gout(13) = gin(15) * gin(18)
      gout(14) = gin(16) * gin(18)
      gout(15) = gin(17) * gin(18)
      gout(16) = gin(19) * gin(24)
      gout(17) = gin(20) * gin(24)
      gout(18) = gin(21) * gin(24)
      gout(19) = gin(22) * gin(24)
      gout(20) = gin(23) * gin(24)
      gout(21) = gin(25) * gin(30)
      gout(22) = gin(26) * gin(30)
      gout(23) = gin(27) * gin(30)
      gout(24) = gin(28) * gin(30)
      gout(25) = gin(29) * gin(30)
      gout(26) = gin(31) * gin(36)
      gout(27) = gin(32) * gin(36)
      gout(28) = gin(33) * gin(36)
      gout(29) = gin(34) * gin(36)
      gout(30) = gin(35) * gin(36)
      gout(31) = gin(37) * gin(42)
      gout(32) = gin(38) * gin(42)
      gout(33) = gin(39) * gin(42)
      gout(34) = gin(40) * gin(42)
      gout(35) = gin(41) * gin(42)
      gout(36) = gin(43) * gin(48)
      gout(37) = gin(44) * gin(48)
      gout(38) = gin(45) * gin(48)
      gout(39) = gin(46) * gin(48)
      gout(40) = gin(47) * gin(48)
    end do
    
    call this%fourtrans%exec_r2c_sub(40, grid, sumNS)
    
  end subroutine grid_op_8_vcst_sub
  
  module pure subroutine grid_op_16_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(96, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+96*i:96+96*i)
      gout => grid(1+80*i:80+80*i)
      
      gout( 1) = gin( 1) * gin( 6)
      gout( 2) = gin( 2) * gin( 6)
      gout( 3) = gin( 3) * gin( 6)
      gout( 4) = gin( 4) * gin( 6)
      gout( 5) = gin( 5) * gin( 6)
      gout( 6) = gin( 7) * gin(12)
      gout( 7) = gin( 8) * gin(12)
      gout( 8) = gin( 9) * gin(12)
      gout( 9) = gin(10) * gin(12)
      gout(10) = gin(11) * gin(12)
      gout(11) = gin(13) * gin(18)
      gout(12) = gin(14) * gin(18)
      gout(13) = gin(15) * gin(18)
      gout(14) = gin(16) * gin(18)
      gout(15) = gin(17) * gin(18)
      gout(16) = gin(19) * gin(24)
      gout(17) = gin(20) * gin(24)
      gout(18) = gin(21) * gin(24)
      gout(19) = gin(22) * gin(24)
      gout(20) = gin(23) * gin(24)
      gout(21) = gin(25) * gin(30)
      gout(22) = gin(26) * gin(30)
      gout(23) = gin(27) * gin(30)
      gout(24) = gin(28) * gin(30)
      gout(25) = gin(29) * gin(30)
      gout(26) = gin(31) * gin(36)
      gout(27) = gin(32) * gin(36)
      gout(28) = gin(33) * gin(36)
      gout(29) = gin(34) * gin(36)
      gout(30) = gin(35) * gin(36)
      gout(31) = gin(37) * gin(42)
      gout(32) = gin(38) * gin(42)
      gout(33) = gin(39) * gin(42)
      gout(34) = gin(40) * gin(42)
      gout(35) = gin(41) * gin(42)
      gout(36) = gin(43) * gin(48)
      gout(37) = gin(44) * gin(48)
      gout(38) = gin(45) * gin(48)
      gout(39) = gin(46) * gin(48)
      gout(40) = gin(47) * gin(48)
      gout(41) = gin(49) * gin(54)
      gout(42) = gin(50) * gin(54)
      gout(43) = gin(51) * gin(54)
      gout(44) = gin(52) * gin(54)
      gout(45) = gin(53) * gin(54)
      gout(46) = gin(55) * gin(60)
      gout(47) = gin(56) * gin(60)
      gout(48) = gin(57) * gin(60)
      gout(49) = gin(58) * gin(60)
      gout(50) = gin(59) * gin(60)
      gout(51) = gin(61) * gin(66)
      gout(52) = gin(62) * gin(66)
      gout(53) = gin(63) * gin(66)
      gout(54) = gin(64) * gin(66)
      gout(55) = gin(65) * gin(66)
      gout(56) = gin(67) * gin(72)
      gout(57) = gin(68) * gin(72)
      gout(58) = gin(69) * gin(72)
      gout(59) = gin(70) * gin(72)
      gout(60) = gin(71) * gin(72)
      gout(61) = gin(73) * gin(78)
      gout(62) = gin(74) * gin(78)
      gout(63) = gin(75) * gin(78)
      gout(64) = gin(76) * gin(78)
      gout(65) = gin(77) * gin(78)
      gout(66) = gin(79) * gin(84)
      gout(67) = gin(80) * gin(84)
      gout(68) = gin(81) * gin(84)
      gout(69) = gin(82) * gin(84)
      gout(70) = gin(83) * gin(84)
      gout(71) = gin(85) * gin(90)
      gout(72) = gin(86) * gin(90)
      gout(73) = gin(87) * gin(90)
      gout(74) = gin(88) * gin(90)
      gout(75) = gin(89) * gin(90)
      gout(76) = gin(91) * gin(96)
      gout(77) = gin(92) * gin(96)
      gout(78) = gin(93) * gin(96)
      gout(79) = gin(94) * gin(96)
      gout(80) = gin(95) * gin(96)
    end do
    
    call this%fourtrans%exec_r2c_sub(80, grid, sumNS)
    
  end subroutine grid_op_16_vcst_sub
  
end submodule vcst_grid