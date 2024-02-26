submodule(SphericalHarmonics) vcvv_vcvgv_grid
  implicit none; contains
  
  module pure subroutine grid_op_2_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(30, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+30*i:30+30*i)
      gout => grid(1+ 8*i: 8+ 8*i)
      
      gout(1) = gin( 1) * gin( 4) + gin( 2) * gin( 5) + gin( 3) * gin( 6)
      gout(2) = gin( 4) * gin( 7) + gin( 5) * gin( 8) + gin( 6) * gin( 9)
      gout(3) = gin( 4) * gin(10) + gin( 5) * gin(11) + gin( 6) * gin(12)
      gout(4) = gin( 4) * gin(13) + gin( 5) * gin(14) + gin( 6) * gin(15)
      gout(5) = gin(16) * gin(19) + gin(17) * gin(20) + gin(18) * gin(21)
      gout(6) = gin(19) * gin(22) + gin(20) * gin(23) + gin(21) * gin(24)
      gout(7) = gin(19) * gin(25) + gin(20) * gin(26) + gin(21) * gin(27)
      gout(8) = gin(19) * gin(28) + gin(20) * gin(29) + gin(21) * gin(30)
    end do
    
    call this%fourtrans%exec_r2c_sub(8, grid, sumNS)
    
  end subroutine grid_op_2_vcvv_vcvgv_sub
  
  module pure subroutine grid_op_4_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(60, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+60*i:60+60*i)
      gout => grid(1+16*i:16+16*i)
      
      gout( 1) = gin( 1) * gin( 4) + gin( 2) * gin( 5) + gin( 3) * gin( 6)
      gout( 2) = gin( 4) * gin( 7) + gin( 5) * gin( 8) + gin( 6) * gin( 9)
      gout( 3) = gin( 4) * gin(10) + gin( 5) * gin(11) + gin( 6) * gin(12)
      gout( 4) = gin( 4) * gin(13) + gin( 5) * gin(14) + gin( 6) * gin(15)
      gout( 5) = gin(16) * gin(19) + gin(17) * gin(20) + gin(18) * gin(21)
      gout( 6) = gin(19) * gin(22) + gin(20) * gin(23) + gin(21) * gin(24)
      gout( 7) = gin(19) * gin(25) + gin(20) * gin(26) + gin(21) * gin(27)
      gout( 8) = gin(19) * gin(28) + gin(20) * gin(29) + gin(21) * gin(30)
      gout( 9) = gin(31) * gin(34) + gin(32) * gin(35) + gin(33) * gin(36)
      gout(10) = gin(34) * gin(37) + gin(35) * gin(38) + gin(36) * gin(39)
      gout(11) = gin(34) * gin(40) + gin(35) * gin(41) + gin(36) * gin(42)
      gout(12) = gin(34) * gin(43) + gin(35) * gin(44) + gin(36) * gin(45)
      gout(13) = gin(46) * gin(49) + gin(47) * gin(50) + gin(48) * gin(51)
      gout(14) = gin(49) * gin(52) + gin(50) * gin(53) + gin(51) * gin(54)
      gout(15) = gin(49) * gin(55) + gin(50) * gin(56) + gin(51) * gin(57)
      gout(16) = gin(49) * gin(58) + gin(50) * gin(59) + gin(51) * gin(60)
    end do
    
    call this%fourtrans%exec_r2c_sub(16, grid, sumNS)
    
  end subroutine grid_op_4_vcvv_vcvgv_sub
  
  module pure subroutine grid_op_8_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(120, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+120*i:120+120*i)
      gout => grid(1+ 32*i: 32+ 32*i)
      
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