submodule (SphericalHarmonics) vcst_grid_4
  implicit none; contains
  
  module pure subroutine grid_op_4_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(24, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+24*i:24+24*i)
      gout => grid(1+20*i:20+20*i)
      
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
    end do
    
    call this%fourtrans%exec_r2c_sub(20, grid, sumNS)
    
  end subroutine grid_op_4_vcst_sub
  
end submodule vcst_grid_4