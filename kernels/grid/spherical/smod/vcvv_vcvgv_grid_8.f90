submodule (SphericalHarmonics) vcvv_vcvgv_grid_8
  implicit none; contains
  
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
  
end submodule vcvv_vcvgv_grid_8