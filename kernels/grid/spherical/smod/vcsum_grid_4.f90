submodule(SphericalHarmonics) vcsum_grid_4
  implicit none; contains
  
  module pure subroutine grid_op_4_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(8, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+8*i:8+8*i)
      gout => grid(1+4*i:4+4*i)
      
      gout(1) = gin(1) * gin(2)
      gout(2) = gin(3) * gin(4)
      gout(3) = gin(5) * gin(6)
      gout(4) = gin(7) * gin(8)
    end do
    
    call this%fourtrans%exec_r2c_sub(4, grid, sumNS)
    
  end subroutine grid_op_4_vcsum_sub
  
end submodule vcsum_grid_4