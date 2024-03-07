submodule (SphericalHarmonics) vcvgv_grid_2
  implicit none; contains
  
  module pure subroutine grid_op_2_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(24, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+24*i:24+24*i)
      gout => grid(1+ 6*i: 6+ 6*i)
      
      gout(1) = gin( 1) * gin(10) + gin( 2) * gin(11) + gin( 3) * gin(12)
      gout(2) = gin( 4) * gin(10) + gin( 5) * gin(11) + gin( 6) * gin(12)
      gout(3) = gin( 7) * gin(10) + gin( 8) * gin(11) + gin( 9) * gin(12)
      gout(4) = gin(13) * gin(22) + gin(14) * gin(23) + gin(15) * gin(24)
      gout(5) = gin(16) * gin(22) + gin(17) * gin(23) + gin(18) * gin(24)
      gout(6) = gin(19) * gin(22) + gin(20) * gin(23) + gin(21) * gin(24)
    end do
    
    call this%fourtrans%exec_r2c_sub(6, grid, sumNS)
    
  end subroutine grid_op_2_vcvgv_sub
  
end submodule vcvgv_grid_2