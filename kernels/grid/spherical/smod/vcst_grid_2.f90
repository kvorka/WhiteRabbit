submodule (SphericalHarmonics) vcst_grid_2
  implicit none; contains
  
  module pure subroutine grid_op_2_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(12, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+12*i:12+12*i)
      gout => grid(1+10*i:10+10*i)
      
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
    end do
    
    call this%fourtrans%exec_r2c_sub(10, grid, sumNS)
    
  end subroutine grid_op_2_vcst_sub
  
end submodule vcst_grid_2