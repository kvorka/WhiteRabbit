submodule (SphericalHarmonics) vcvv_vcvgv_grid_4
  implicit none; contains
  
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
  
end submodule vcvv_vcvgv_grid_4