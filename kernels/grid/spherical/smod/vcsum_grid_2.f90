submodule(SphericalHarmonics) vcsum_grid_2
  implicit none; contains
  
  module pure subroutine grid_op_2_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    
    call this%fourtrans%exec_c2r_sub(4, sumNS, grid)
    
    do i = 0, this%nFourier-1
      gin  => grid(1+4*i:4+4*i)
      gout => grid(1+2*i:2+2*i)
      
      gout(1) = gin(1) * gin(2)
      gout(2) = gin(3) * gin(4)
    end do
    
    call this%fourtrans%exec_r2c_sub(2, grid, sumNS)
    
  end subroutine grid_op_2_vcsum_sub
  
end submodule vcsum_grid_2