submodule (SphericalHarmonics) vcvv_vcvgv_grid_16
  implicit none; contains
  
  module pure subroutine grid_op_16_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gout(:), gin(:)
    real(kind=dbl), allocatable           :: tmp(:)
    
    call this%fourtrans%exec_c2r_sub(240, sumNS, grid)
    
    allocate( tmp(3) )
    
    do i = 0, this%nFourier-1
      gin  => grid(1+240*i:240+240*i)
      gout => grid(1+ 64*i: 64+ 64*i)
      
      tmp = gin(1:3)
        gout(1) = sum( gin( 4: 6) * tmp )
        gout(2) = sum( gin( 7: 9) * tmp )
        gout(3) = sum( gin(10:12) * tmp )
        gout(4) = sum( gin(13:15) * tmp )
      
      tmp = gin(16:18)
        gout(5) = sum( gin(19:21) * tmp )
        gout(6) = sum( gin(22:24) * tmp )
        gout(7) = sum( gin(25:27) * tmp )
        gout(8) = sum( gin(28:30) * tmp )
      
      tmp = gin(31:33)
        gout( 9) = sum( gin(34:36) * tmp )
        gout(10) = sum( gin(37:39) * tmp )
        gout(11) = sum( gin(40:42) * tmp )
        gout(12) = sum( gin(43:45) * tmp )
      
      tmp = gin(46:48)
        gout(13) = sum( gin(49:51) * tmp )
        gout(14) = sum( gin(52:54) * tmp )
        gout(15) = sum( gin(55:57) * tmp )
        gout(16) = sum( gin(58:60) * tmp )
      
      tmp = gin(61:63)
        gout(17) = sum( gin(64:66) * tmp )
        gout(18) = sum( gin(67:69) * tmp )
        gout(19) = sum( gin(70:72) * tmp )
        gout(20) = sum( gin(73:75) * tmp )
      
      tmp = gin(76:78)
        gout(21) = sum( gin(79:81) * tmp )
        gout(22) = sum( gin(82:84) * tmp )
        gout(23) = sum( gin(85:87) * tmp )
        gout(24) = sum( gin(88:90) * tmp )
      
      tmp = gin(91:93)
        gout(25) = sum( gin( 94: 96) * tmp )
        gout(26) = sum( gin( 97: 99) * tmp )
        gout(27) = sum( gin(100:102) * tmp )
        gout(28) = sum( gin(103:105) * tmp )
      
      tmp = gin(106:108)
        gout(29) = sum( gin(109:111) * tmp )
        gout(30) = sum( gin(112:114) * tmp )
        gout(31) = sum( gin(115:117) * tmp )
        gout(32) = sum( gin(118:120) * tmp )
      
      tmp = gin(121:123)
        gout(33) = sum( gin(124:126) * tmp )
        gout(34) = sum( gin(127:129) * tmp )
        gout(35) = sum( gin(130:132) * tmp )
        gout(36) = sum( gin(133:135) * tmp )
      
      tmp = gin(136:138)
        gout(37) = sum( gin(139:141) * tmp )
        gout(38) = sum( gin(142:144) * tmp )
        gout(39) = sum( gin(145:147) * tmp )
        gout(40) = sum( gin(148:150) * tmp )
      
      tmp = gin(151:153)
        gout(41) = sum( gin(154:156) * tmp )
        gout(42) = sum( gin(157:159) * tmp )
        gout(43) = sum( gin(160:162) * tmp )
        gout(44) = sum( gin(163:165) * tmp )
      
      tmp = gin(166:168)
        gout(45) = sum( gin(169:171) * tmp )
        gout(46) = sum( gin(172:174) * tmp )
        gout(47) = sum( gin(175:177) * tmp )
        gout(48) = sum( gin(178:180) * tmp )
      
      tmp = gin(181:183)
        gout(49) = sum( gin(184:186) * tmp )
        gout(50) = sum( gin(187:189) * tmp )
        gout(51) = sum( gin(190:192) * tmp )
        gout(52) = sum( gin(192:195) * tmp )
      
      tmp = gin(196:198)
        gout(53) = sum( gin(199:201) * tmp )
        gout(54) = sum( gin(202:204) * tmp )
        gout(55) = sum( gin(205:207) * tmp )
        gout(56) = sum( gin(208:210) * tmp )
      
      tmp = gin(211:213)
        gout(57) = sum( gin(214:216) * tmp )
        gout(58) = sum( gin(217:219) * tmp )
        gout(59) = sum( gin(220:222) * tmp )
        gout(60) = sum( gin(223:225) * tmp )
      
      tmp = gin(226:228)
        gout(61) = sum( gin(229:231) * tmp )
        gout(62) = sum( gin(232:234) * tmp )
        gout(63) = sum( gin(235:237) * tmp )
        gout(64) = sum( gin(238:240) * tmp )
    end do
    
    deallocate( tmp )
    
    call this%fourtrans%exec_r2c_sub(64, grid, sumNS)
    
  end subroutine grid_op_16_vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv_grid_16