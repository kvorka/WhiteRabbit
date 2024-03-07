submodule (SphericalHarmonics) vcvgv_grid_8
  implicit none; contains
  
  module pure subroutine grid_op_8_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(96, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+96*i:96+96*i)
      gout => grid(1+24*i:24+24*i)
      
      gout( 1) = gin( 1) * gin(10) + gin( 2) * gin(11) + gin( 3) * gin(12)
      gout( 2) = gin( 4) * gin(10) + gin( 5) * gin(11) + gin( 6) * gin(12)
      gout( 3) = gin( 7) * gin(10) + gin( 8) * gin(11) + gin( 9) * gin(12)
      gout( 4) = gin(13) * gin(22) + gin(14) * gin(23) + gin(15) * gin(24)
      gout( 5) = gin(16) * gin(22) + gin(17) * gin(23) + gin(18) * gin(24)
      gout( 6) = gin(19) * gin(22) + gin(20) * gin(23) + gin(21) * gin(24)
      gout( 7) = gin(25) * gin(34) + gin(26) * gin(35) + gin(27) * gin(36)
      gout( 8) = gin(28) * gin(34) + gin(29) * gin(35) + gin(30) * gin(36)
      gout( 9) = gin(31) * gin(34) + gin(32) * gin(35) + gin(33) * gin(36)
      gout(10) = gin(37) * gin(46) + gin(38) * gin(47) + gin(39) * gin(48)
      gout(11) = gin(40) * gin(46) + gin(41) * gin(47) + gin(42) * gin(48)
      gout(12) = gin(43) * gin(46) + gin(44) * gin(47) + gin(45) * gin(48)
      gout(13) = gin(49) * gin(58) + gin(50) * gin(59) + gin(51) * gin(60)
      gout(14) = gin(52) * gin(58) + gin(53) * gin(59) + gin(54) * gin(60)
      gout(15) = gin(55) * gin(58) + gin(56) * gin(59) + gin(57) * gin(60)
      gout(16) = gin(61) * gin(70) + gin(62) * gin(71) + gin(63) * gin(72)
      gout(17) = gin(64) * gin(70) + gin(65) * gin(71) + gin(66) * gin(72)
      gout(18) = gin(67) * gin(70) + gin(68) * gin(71) + gin(69) * gin(72)
      gout(19) = gin(73) * gin(82) + gin(74) * gin(83) + gin(75) * gin(84)
      gout(20) = gin(76) * gin(82) + gin(77) * gin(83) + gin(78) * gin(84)
      gout(21) = gin(79) * gin(82) + gin(80) * gin(83) + gin(81) * gin(84)
      gout(22) = gin(85) * gin(94) + gin(86) * gin(95) + gin(87) * gin(96)
      gout(23) = gin(88) * gin(94) + gin(89) * gin(95) + gin(90) * gin(96)
      gout(24) = gin(91) * gin(94) + gin(92) * gin(95) + gin(93) * gin(96)
    end do
    
    call this%fourtrans%exec_r2c_sub(24, grid, sumNS)
    
  end subroutine grid_op_8_vcvgv_sub
  
end submodule vcvgv_grid_8