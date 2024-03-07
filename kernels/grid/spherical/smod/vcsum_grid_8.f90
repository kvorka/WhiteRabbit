submodule(SphericalHarmonics) vcsum_grid_8
  implicit none; contains
  
  module pure subroutine grid_op_8_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(16, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+16*i:16+16*i)
      gout => grid(1+ 8*i: 8+ 8*i)
      
      gout(1) = gin( 1) * gin( 2)
      gout(2) = gin( 3) * gin( 4)
      gout(3) = gin( 5) * gin( 6)
      gout(4) = gin( 7) * gin( 8)
      gout(5) = gin( 9) * gin(10)
      gout(6) = gin(11) * gin(12)
      gout(7) = gin(13) * gin(14)
      gout(8) = gin(15) * gin(16)
    end do
    
    call this%fourtrans%exec_r2c_sub(8, grid, sumNS)
    
  end subroutine grid_op_8_vcsum_sub
  
end submodule vcsum_grid_8