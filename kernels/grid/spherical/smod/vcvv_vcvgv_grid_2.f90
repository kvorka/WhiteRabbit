submodule (SphericalHarmonics) vcvv_vcvgv_grid_2
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
  
end submodule vcvv_vcvgv_grid_2