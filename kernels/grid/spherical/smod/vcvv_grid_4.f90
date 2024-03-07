submodule (SphericalHarmonics) vcvv_grid_4
  implicit none; contains
  
  module pure subroutine grid_op_4_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(24, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+24*i:24+24*i)
      gout => grid(1+ 4*i: 4+ 4*i)
      
      gout(1) = gin( 1) * gin( 4) + gin( 2) * gin( 5) + gin( 3) * gin( 6)
      gout(2) = gin( 7) * gin(10) + gin( 8) * gin(11) + gin( 9) * gin(12)
      gout(3) = gin(13) * gin(16) + gin(14) * gin(17) + gin(15) * gin(18)
      gout(4) = gin(19) * gin(22) + gin(20) * gin(23) + gin(21) * gin(24)
    end do
    
    call this%fourtrans%exec_r2c_sub(4, grid, sumNS)
    
  end subroutine grid_op_4_vcvv_sub
  
end submodule vcvv_grid_4