submodule(SphericalHarmonics) vcvv_vcvgv_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i
    
    call this%fourtrans%exec_c2r_sub(30, sumNS, grid)
    
    do i = 0, this%nFourier-1
      grid(1+8*i) = grid( 1+30*i) * grid( 4+30*i) + grid( 2+30*i) * grid( 5+30*i) + grid( 3+30*i) * grid( 6+30*i)
      grid(2+8*i) = grid( 4+30*i) * grid( 7+30*i) + grid( 5+30*i) * grid( 8+30*i) + grid( 6+30*i) * grid( 9+30*i)
      grid(3+8*i) = grid( 4+30*i) * grid(10+30*i) + grid( 5+30*i) * grid(11+30*i) + grid( 6+30*i) * grid(12+30*i)
      grid(4+8*i) = grid( 4+30*i) * grid(13+30*i) + grid( 5+30*i) * grid(14+30*i) + grid( 6+30*i) * grid(15+30*i)
      grid(5+8*i) = grid(16+30*i) * grid(19+30*i) + grid(17+30*i) * grid(20+30*i) + grid(18+30*i) * grid(21+30*i)
      grid(6+8*i) = grid(19+30*i) * grid(22+30*i) + grid(20+30*i) * grid(23+30*i) + grid(21+30*i) * grid(24+30*i)
      grid(7+8*i) = grid(19+30*i) * grid(25+30*i) + grid(20+30*i) * grid(26+30*i) + grid(21+30*i) * grid(27+30*i)
      grid(8+8*i) = grid(19+30*i) * grid(28+30*i) + grid(20+30*i) * grid(29+30*i) + grid(21+30*i) * grid(30+30*i)
    end do
    
    call this%fourtrans%exec_r2c_sub(8, grid, sumNS)
    
  end subroutine grid_op_2_vcvv_vcvgv_sub
  
  module pure subroutine grid_op_4_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i
    
    call this%fourtrans%exec_c2r_sub(60, sumNS, grid)
    
    do i = 0, this%nFourier-1
      grid( 1+16*i) = grid( 1+60*i) * grid( 4+60*i) + grid( 2+60*i) * grid( 5+60*i) + grid( 3+60*i) * grid( 6+60*i)
      grid( 2+16*i) = grid( 4+60*i) * grid( 7+60*i) + grid( 5+60*i) * grid( 8+60*i) + grid( 6+60*i) * grid( 9+60*i)
      grid( 3+16*i) = grid( 4+60*i) * grid(10+60*i) + grid( 5+60*i) * grid(11+60*i) + grid( 6+60*i) * grid(12+60*i)
      grid( 4+16*i) = grid( 4+60*i) * grid(13+60*i) + grid( 5+60*i) * grid(14+60*i) + grid( 6+60*i) * grid(15+60*i)
      grid( 5+16*i) = grid(16+60*i) * grid(19+60*i) + grid(17+60*i) * grid(20+60*i) + grid(18+60*i) * grid(21+60*i)
      grid( 6+16*i) = grid(19+60*i) * grid(22+60*i) + grid(20+60*i) * grid(23+60*i) + grid(21+60*i) * grid(24+60*i)
      grid( 7+16*i) = grid(19+60*i) * grid(25+60*i) + grid(20+60*i) * grid(26+60*i) + grid(21+60*i) * grid(27+60*i)
      grid( 8+16*i) = grid(19+60*i) * grid(28+60*i) + grid(20+60*i) * grid(29+60*i) + grid(21+60*i) * grid(30+60*i)
      grid( 9+16*i) = grid(31+60*i) * grid(34+60*i) + grid(32+60*i) * grid(35+60*i) + grid(33+60*i) * grid(36+60*i)
      grid(10+16*i) = grid(34+60*i) * grid(37+60*i) + grid(35+60*i) * grid(38+60*i) + grid(36+60*i) * grid(39+60*i)
      grid(11+16*i) = grid(34+60*i) * grid(40+60*i) + grid(35+60*i) * grid(41+60*i) + grid(36+60*i) * grid(42+60*i)
      grid(12+16*i) = grid(34+60*i) * grid(43+60*i) + grid(35+60*i) * grid(44+60*i) + grid(36+60*i) * grid(45+60*i)
      grid(13+16*i) = grid(46+60*i) * grid(49+60*i) + grid(47+60*i) * grid(50+60*i) + grid(48+60*i) * grid(51+60*i)
      grid(14+16*i) = grid(49+60*i) * grid(52+60*i) + grid(50+60*i) * grid(53+60*i) + grid(51+60*i) * grid(54+60*i)
      grid(15+16*i) = grid(49+60*i) * grid(55+60*i) + grid(50+60*i) * grid(56+60*i) + grid(51+60*i) * grid(57+60*i)
      grid(16+16*i) = grid(49+60*i) * grid(58+60*i) + grid(50+60*i) * grid(59+60*i) + grid(51+60*i) * grid(60+60*i)
    end do
    
    call this%fourtrans%exec_r2c_sub(16, grid, sumNS)
    
  end subroutine grid_op_4_vcvv_vcvgv_sub
  
  module pure subroutine grid_op_8_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(out)   :: grid(*)
    complex(kind=dbl),    intent(inout) :: sumNS(*)
    integer                             :: i
    
    call this%fourtrans%exec_c2r_sub(120, sumNS, grid)
    
    do i = 0, this%nFourier-1
      grid( 1+32*i) = grid(  1+120*i) * grid(  4+120*i) + grid(  2+120*i) * grid(  5+120*i) + grid(  3+120*i) * grid(  6+120*i)
      grid( 2+32*i) = grid(  4+120*i) * grid(  7+120*i) + grid(  5+120*i) * grid(  8+120*i) + grid(  6+120*i) * grid(  9+120*i)
      grid( 3+32*i) = grid(  4+120*i) * grid( 10+120*i) + grid(  5+120*i) * grid( 11+120*i) + grid(  6+120*i) * grid( 12+120*i)
      grid( 4+32*i) = grid(  4+120*i) * grid( 13+120*i) + grid(  5+120*i) * grid( 14+120*i) + grid(  6+120*i) * grid( 15+120*i)
      grid( 5+32*i) = grid( 16+120*i) * grid( 19+120*i) + grid( 17+120*i) * grid( 20+120*i) + grid( 18+120*i) * grid( 21+120*i)
      grid( 6+32*i) = grid( 19+120*i) * grid( 22+120*i) + grid( 20+120*i) * grid( 23+120*i) + grid( 21+120*i) * grid( 24+120*i)
      grid( 7+32*i) = grid( 19+120*i) * grid( 25+120*i) + grid( 20+120*i) * grid( 26+120*i) + grid( 21+120*i) * grid( 27+120*i)
      grid( 8+32*i) = grid( 19+120*i) * grid( 28+120*i) + grid( 20+120*i) * grid( 29+120*i) + grid( 21+120*i) * grid( 30+120*i)
      grid( 9+32*i) = grid( 31+120*i) * grid( 34+120*i) + grid( 32+120*i) * grid( 35+120*i) + grid( 33+120*i) * grid( 36+120*i)
      grid(10+32*i) = grid( 34+120*i) * grid( 37+120*i) + grid( 35+120*i) * grid( 38+120*i) + grid( 36+120*i) * grid( 39+120*i)
      grid(11+32*i) = grid( 34+120*i) * grid( 40+120*i) + grid( 35+120*i) * grid( 41+120*i) + grid( 36+120*i) * grid( 42+120*i)
      grid(12+32*i) = grid( 34+120*i) * grid( 43+120*i) + grid( 35+120*i) * grid( 44+120*i) + grid( 36+120*i) * grid( 45+120*i)
      grid(13+32*i) = grid( 46+120*i) * grid( 49+120*i) + grid( 47+120*i) * grid( 50+120*i) + grid( 48+120*i) * grid( 51+120*i)
      grid(14+32*i) = grid( 49+120*i) * grid( 52+120*i) + grid( 50+120*i) * grid( 53+120*i) + grid( 51+120*i) * grid( 54+120*i)
      grid(15+32*i) = grid( 49+120*i) * grid( 55+120*i) + grid( 50+120*i) * grid( 56+120*i) + grid( 51+120*i) * grid( 57+120*i)
      grid(16+32*i) = grid( 49+120*i) * grid( 58+120*i) + grid( 50+120*i) * grid( 59+120*i) + grid( 51+120*i) * grid( 60+120*i)
      grid(17+32*i) = grid( 61+120*i) * grid( 64+120*i) + grid( 62+120*i) * grid( 65+120*i) + grid( 63+120*i) * grid( 66+120*i)
      grid(18+32*i) = grid( 64+120*i) * grid( 67+120*i) + grid( 65+120*i) * grid( 68+120*i) + grid( 66+120*i) * grid( 69+120*i)
      grid(19+32*i) = grid( 64+120*i) * grid( 70+120*i) + grid( 65+120*i) * grid( 71+120*i) + grid( 66+120*i) * grid( 72+120*i)
      grid(20+32*i) = grid( 64+120*i) * grid( 73+120*i) + grid( 65+120*i) * grid( 74+120*i) + grid( 66+120*i) * grid( 75+120*i)
      grid(21+32*i) = grid( 76+120*i) * grid( 79+120*i) + grid( 77+120*i) * grid( 80+120*i) + grid( 78+120*i) * grid( 81+120*i)
      grid(22+32*i) = grid( 79+120*i) * grid( 82+120*i) + grid( 80+120*i) * grid( 83+120*i) + grid( 81+120*i) * grid( 84+120*i)
      grid(23+32*i) = grid( 79+120*i) * grid( 85+120*i) + grid( 80+120*i) * grid( 86+120*i) + grid( 81+120*i) * grid( 87+120*i)
      grid(24+32*i) = grid( 79+120*i) * grid( 88+120*i) + grid( 80+120*i) * grid( 89+120*i) + grid( 81+120*i) * grid( 90+120*i)
      grid(25+32*i) = grid( 91+120*i) * grid( 94+120*i) + grid( 92+120*i) * grid( 95+120*i) + grid( 93+120*i) * grid( 96+120*i)
      grid(26+32*i) = grid( 94+120*i) * grid( 97+120*i) + grid( 95+120*i) * grid( 98+120*i) + grid( 96+120*i) * grid( 99+120*i)
      grid(27+32*i) = grid( 94+120*i) * grid(100+120*i) + grid( 95+120*i) * grid(101+120*i) + grid( 96+120*i) * grid(102+120*i)
      grid(28+32*i) = grid( 94+120*i) * grid(103+120*i) + grid( 95+120*i) * grid(104+120*i) + grid( 96+120*i) * grid(105+120*i)
      grid(29+32*i) = grid(106+120*i) * grid(109+120*i) + grid(107+120*i) * grid(110+120*i) + grid(108+120*i) * grid(111+120*i)
      grid(30+32*i) = grid(109+120*i) * grid(112+120*i) + grid(110+120*i) * grid(113+120*i) + grid(111+120*i) * grid(114+120*i)
      grid(31+32*i) = grid(109+120*i) * grid(115+120*i) + grid(110+120*i) * grid(116+120*i) + grid(111+120*i) * grid(117+120*i)
      grid(32+32*i) = grid(109+120*i) * grid(118+120*i) + grid(110+120*i) * grid(119+120*i) + grid(111+120*i) * grid(120+120*i)
    end do
    
    call this%fourtrans%exec_r2c_sub(32, grid, sumNS)
    
  end subroutine grid_op_8_vcvv_vcvgv_sub
  
  module pure subroutine grid_op_16_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gout(:), gin(:)
    
    call this%fourtrans%exec_c2r_sub(240, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+240*i:240+240*i)
      gout => grid(1+ 64*i: 64+ 64*i)
      
      gout( 1) = gin(  1) * gin(  4) + gin(  2) * gin(  5) + gin(  3) * gin(  6)
      gout( 2) = gin(  4) * gin(  7) + gin(  5) * gin(  8) + gin(  6) * gin(  9)
      gout( 3) = gin(  4) * gin( 10) + gin(  5) * gin( 11) + gin(  6) * gin( 12)
      gout( 4) = gin(  4) * gin( 13) + gin(  5) * gin( 14) + gin(  6) * gin( 15)
      gout( 5) = gin( 16) * gin( 19) + gin( 17) * gin( 20) + gin( 18) * gin( 21)
      gout( 6) = gin( 19) * gin( 22) + gin( 20) * gin( 23) + gin( 21) * gin( 24)
      gout( 7) = gin( 19) * gin( 25) + gin( 20) * gin( 26) + gin( 21) * gin( 27)
      gout( 8) = gin( 19) * gin( 28) + gin( 20) * gin( 29) + gin( 21) * gin( 30)
      gout( 9) = gin( 31) * gin( 34) + gin( 32) * gin( 35) + gin( 33) * gin( 36)
      gout(10) = gin( 34) * gin( 37) + gin( 35) * gin( 38) + gin( 36) * gin( 39)
      gout(11) = gin( 34) * gin( 40) + gin( 35) * gin( 41) + gin( 36) * gin( 42)
      gout(12) = gin( 34) * gin( 43) + gin( 35) * gin( 44) + gin( 36) * gin( 45)
      gout(13) = gin( 46) * gin( 49) + gin( 47) * gin( 50) + gin( 48) * gin( 51)
      gout(14) = gin( 49) * gin( 52) + gin( 50) * gin( 53) + gin( 51) * gin( 54)
      gout(15) = gin( 49) * gin( 55) + gin( 50) * gin( 56) + gin( 51) * gin( 57)
      gout(16) = gin( 49) * gin( 58) + gin( 50) * gin( 59) + gin( 51) * gin( 60)
      gout(17) = gin( 61) * gin( 64) + gin( 62) * gin( 65) + gin( 63) * gin( 66)
      gout(18) = gin( 64) * gin( 67) + gin( 65) * gin( 68) + gin( 66) * gin( 69)
      gout(19) = gin( 64) * gin( 70) + gin( 65) * gin( 71) + gin( 66) * gin( 72)
      gout(20) = gin( 64) * gin( 73) + gin( 65) * gin( 74) + gin( 66) * gin( 75)
      gout(21) = gin( 76) * gin( 79) + gin( 77) * gin( 80) + gin( 78) * gin( 81)
      gout(22) = gin( 79) * gin( 82) + gin( 80) * gin( 83) + gin( 81) * gin( 84)
      gout(23) = gin( 79) * gin( 85) + gin( 80) * gin( 86) + gin( 81) * gin( 87)
      gout(24) = gin( 79) * gin( 88) + gin( 80) * gin( 89) + gin( 81) * gin( 90)
      gout(25) = gin( 91) * gin( 94) + gin( 92) * gin( 95) + gin( 93) * gin( 96)
      gout(26) = gin( 94) * gin( 97) + gin( 95) * gin( 98) + gin( 96) * gin( 99)
      gout(27) = gin( 94) * gin(100) + gin( 95) * gin(101) + gin( 96) * gin(102)
      gout(28) = gin( 94) * gin(103) + gin( 95) * gin(104) + gin( 96) * gin(105)
      gout(29) = gin(106) * gin(109) + gin(107) * gin(110) + gin(108) * gin(111)
      gout(30) = gin(109) * gin(112) + gin(110) * gin(113) + gin(111) * gin(114)
      gout(31) = gin(109) * gin(115) + gin(110) * gin(116) + gin(111) * gin(117)
      gout(32) = gin(109) * gin(118) + gin(110) * gin(119) + gin(111) * gin(120)
      gout(33) = gin(121) * gin(124) + gin(122) * gin(125) + gin(123) * gin(126)
      gout(34) = gin(124) * gin(127) + gin(125) * gin(128) + gin(126) * gin(129)
      gout(35) = gin(124) * gin(130) + gin(125) * gin(131) + gin(126) * gin(132)
      gout(36) = gin(124) * gin(133) + gin(125) * gin(134) + gin(126) * gin(135)
      gout(37) = gin(136) * gin(139) + gin(137) * gin(140) + gin(138) * gin(141)
      gout(38) = gin(139) * gin(142) + gin(140) * gin(143) + gin(141) * gin(144)
      gout(39) = gin(139) * gin(145) + gin(140) * gin(146) + gin(141) * gin(147)
      gout(40) = gin(139) * gin(148) + gin(140) * gin(149) + gin(141) * gin(150)
      gout(41) = gin(151) * gin(154) + gin(152) * gin(155) + gin(153) * gin(156)
      gout(42) = gin(154) * gin(157) + gin(155) * gin(158) + gin(156) * gin(159)
      gout(43) = gin(154) * gin(160) + gin(155) * gin(161) + gin(156) * gin(162)
      gout(44) = gin(154) * gin(163) + gin(155) * gin(164) + gin(156) * gin(165)
      gout(45) = gin(166) * gin(169) + gin(167) * gin(170) + gin(168) * gin(171)
      gout(46) = gin(169) * gin(172) + gin(170) * gin(173) + gin(171) * gin(174)
      gout(47) = gin(169) * gin(175) + gin(170) * gin(176) + gin(171) * gin(177)
      gout(48) = gin(169) * gin(178) + gin(170) * gin(179) + gin(171) * gin(180)
      gout(49) = gin(181) * gin(184) + gin(182) * gin(185) + gin(183) * gin(186)
      gout(50) = gin(184) * gin(187) + gin(185) * gin(188) + gin(186) * gin(189)
      gout(51) = gin(184) * gin(190) + gin(185) * gin(191) + gin(186) * gin(192)
      gout(52) = gin(184) * gin(193) + gin(185) * gin(194) + gin(186) * gin(195)
      gout(53) = gin(196) * gin(199) + gin(197) * gin(200) + gin(198) * gin(201)
      gout(54) = gin(199) * gin(202) + gin(200) * gin(203) + gin(201) * gin(204)
      gout(55) = gin(199) * gin(205) + gin(200) * gin(206) + gin(201) * gin(207)
      gout(56) = gin(199) * gin(208) + gin(200) * gin(209) + gin(201) * gin(210)
      gout(57) = gin(211) * gin(214) + gin(212) * gin(215) + gin(213) * gin(216)
      gout(58) = gin(214) * gin(217) + gin(215) * gin(218) + gin(216) * gin(219)
      gout(59) = gin(214) * gin(220) + gin(215) * gin(221) + gin(216) * gin(222)
      gout(60) = gin(214) * gin(223) + gin(215) * gin(224) + gin(216) * gin(225)
      gout(61) = gin(226) * gin(229) + gin(227) * gin(230) + gin(228) * gin(231)
      gout(62) = gin(229) * gin(232) + gin(230) * gin(233) + gin(231) * gin(234)
      gout(63) = gin(229) * gin(235) + gin(230) * gin(236) + gin(231) * gin(237)
      gout(64) = gin(229) * gin(238) + gin(230) * gin(239) + gin(231) * gin(240)
    end do
    
    call this%fourtrans%exec_r2c_sub(64, grid, sumNS)
    
  end subroutine grid_op_16_vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv_grid