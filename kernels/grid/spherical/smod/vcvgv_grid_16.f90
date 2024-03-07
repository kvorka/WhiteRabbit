submodule (SphericalHarmonics) vcvgv_grid_16
  implicit none; contains
  
  module pure subroutine grid_op_16_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(192, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+192*i:192+192*i)
      gout => grid(1+ 48*i: 48+ 48*i)
      
      gout( 1) = gin(  1) * gin( 10) + gin(  2) * gin( 11) + gin(  3) * gin( 12)
      gout( 2) = gin(  4) * gin( 10) + gin(  5) * gin( 11) + gin(  6) * gin( 12)
      gout( 3) = gin(  7) * gin( 10) + gin(  8) * gin( 11) + gin(  9) * gin( 12)
      gout( 4) = gin( 13) * gin( 22) + gin( 14) * gin( 23) + gin( 15) * gin( 24)
      gout( 5) = gin( 16) * gin( 22) + gin( 17) * gin( 23) + gin( 18) * gin( 24)
      gout( 6) = gin( 19) * gin( 22) + gin( 20) * gin( 23) + gin( 21) * gin( 24)
      gout( 7) = gin( 25) * gin( 34) + gin( 26) * gin( 35) + gin( 27) * gin( 36)
      gout( 8) = gin( 28) * gin( 34) + gin( 29) * gin( 35) + gin( 30) * gin( 36)
      gout( 9) = gin( 31) * gin( 34) + gin( 32) * gin( 35) + gin( 33) * gin( 36)
      gout(10) = gin( 37) * gin( 46) + gin( 38) * gin( 47) + gin( 39) * gin( 48)
      gout(11) = gin( 40) * gin( 46) + gin( 41) * gin( 47) + gin( 42) * gin( 48)
      gout(12) = gin( 43) * gin( 46) + gin( 44) * gin( 47) + gin( 45) * gin( 48)
      gout(13) = gin( 49) * gin( 58) + gin( 50) * gin( 59) + gin( 51) * gin( 60)
      gout(14) = gin( 52) * gin( 58) + gin( 53) * gin( 59) + gin( 54) * gin( 60)
      gout(15) = gin( 55) * gin( 58) + gin( 56) * gin( 59) + gin( 57) * gin( 60)
      gout(16) = gin( 61) * gin( 70) + gin( 62) * gin( 71) + gin( 63) * gin( 72)
      gout(17) = gin( 64) * gin( 70) + gin( 65) * gin( 71) + gin( 66) * gin( 72)
      gout(18) = gin( 67) * gin( 70) + gin( 68) * gin( 71) + gin( 69) * gin( 72)
      gout(19) = gin( 73) * gin( 82) + gin( 74) * gin( 83) + gin( 75) * gin( 84)
      gout(20) = gin( 76) * gin( 82) + gin( 77) * gin( 83) + gin( 78) * gin( 84)
      gout(21) = gin( 79) * gin( 82) + gin( 80) * gin( 83) + gin( 81) * gin( 84)
      gout(22) = gin( 85) * gin( 94) + gin( 86) * gin( 95) + gin( 87) * gin( 96)
      gout(23) = gin( 88) * gin( 94) + gin( 89) * gin( 95) + gin( 90) * gin( 96)
      gout(24) = gin( 91) * gin( 94) + gin( 92) * gin( 95) + gin( 93) * gin( 96)
      gout(25) = gin( 97) * gin(106) + gin( 98) * gin(107) + gin( 99) * gin(108)
      gout(26) = gin(100) * gin(106) + gin(101) * gin(107) + gin(102) * gin(108)
      gout(27) = gin(103) * gin(106) + gin(104) * gin(107) + gin(105) * gin(108)
      gout(28) = gin(109) * gin(118) + gin(110) * gin(119) + gin(111) * gin(120)
      gout(29) = gin(112) * gin(118) + gin(113) * gin(119) + gin(114) * gin(120)
      gout(30) = gin(115) * gin(118) + gin(116) * gin(119) + gin(117) * gin(120)
      gout(31) = gin(121) * gin(130) + gin(122) * gin(131) + gin(123) * gin(132)
      gout(32) = gin(124) * gin(130) + gin(125) * gin(131) + gin(126) * gin(132)
      gout(33) = gin(127) * gin(130) + gin(128) * gin(131) + gin(129) * gin(132)
      gout(34) = gin(133) * gin(142) + gin(134) * gin(143) + gin(135) * gin(144)
      gout(35) = gin(136) * gin(142) + gin(137) * gin(143) + gin(138) * gin(144)
      gout(36) = gin(139) * gin(142) + gin(140) * gin(143) + gin(141) * gin(144)
      gout(37) = gin(145) * gin(154) + gin(146) * gin(155) + gin(147) * gin(156)
      gout(38) = gin(148) * gin(154) + gin(149) * gin(155) + gin(150) * gin(156)
      gout(39) = gin(151) * gin(154) + gin(152) * gin(155) + gin(153) * gin(156)
      gout(40) = gin(157) * gin(166) + gin(158) * gin(167) + gin(159) * gin(168)
      gout(41) = gin(160) * gin(166) + gin(161) * gin(167) + gin(162) * gin(168)
      gout(42) = gin(163) * gin(166) + gin(164) * gin(167) + gin(165) * gin(168)
      gout(43) = gin(169) * gin(178) + gin(170) * gin(179) + gin(171) * gin(180)
      gout(44) = gin(172) * gin(178) + gin(173) * gin(179) + gin(174) * gin(180)
      gout(45) = gin(175) * gin(178) + gin(176) * gin(179) + gin(177) * gin(180)
      gout(46) = gin(181) * gin(190) + gin(182) * gin(191) + gin(183) * gin(192)
      gout(47) = gin(184) * gin(190) + gin(185) * gin(191) + gin(186) * gin(192)
      gout(48) = gin(187) * gin(190) + gin(188) * gin(191) + gin(189) * gin(192)
    end do
    
    call this%fourtrans%exec_r2c_sub(48, grid, sumNS)
    
  end subroutine grid_op_16_vcvgv_sub
  
end submodule vcvgv_grid_16