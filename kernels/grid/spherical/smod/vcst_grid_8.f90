submodule (SphericalHarmonics) vcst_grid_8
  implicit none; contains
  
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
  
end submodule vcst_grid_8