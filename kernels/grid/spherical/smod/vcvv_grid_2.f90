submodule (SphericalHarmonics) vcvv_grid_2
  implicit none; contains
  
  module pure subroutine grid_op_2_vcvv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(12, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+12*i:12+12*i)
      gout => grid(1+ 2*i: 2+ 2*i)
      
      gout(1) = gin(1) * gin( 4) + gin(2) * gin( 5) + gin(3) * gin( 6)
      gout(2) = gin(7) * gin(10) + gin(8) * gin(11) + gin(9) * gin(12)
    end do
    
    call this%fourtrans%exec_r2c_sub(2, grid, sumNS)
    
  end subroutine grid_op_2_vcvv_sub
  
end submodule vcvv_grid_2