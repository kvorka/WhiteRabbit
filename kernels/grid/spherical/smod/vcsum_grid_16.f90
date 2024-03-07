submodule(SphericalHarmonics) vcsum_grid_16
  implicit none; contains
  
  module pure subroutine grid_op_16_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(32, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+32*i:32+32*i)
      gout => grid(1+16*i:16+16*i)
      
      gout( 1) = gin( 1) * gin( 2)
      gout( 2) = gin( 3) * gin( 4)
      gout( 3) = gin( 5) * gin( 6)
      gout( 4) = gin( 7) * gin( 8)
      gout( 5) = gin( 9) * gin(10)
      gout( 6) = gin(11) * gin(12)
      gout( 7) = gin(13) * gin(14)
      gout( 8) = gin(15) * gin(16)
      gout( 9) = gin(17) * gin(18)
      gout(10) = gin(19) * gin(20)
      gout(11) = gin(21) * gin(22)
      gout(12) = gin(23) * gin(24)
      gout(13) = gin(25) * gin(26)
      gout(14) = gin(27) * gin(28)
      gout(15) = gin(29) * gin(30)
      gout(16) = gin(31) * gin(32)
    end do
    
    call this%fourtrans%exec_r2c_sub(16, grid, sumNS)
    
  end subroutine grid_op_16_vcsum_sub
  
end submodule vcsum_grid_16