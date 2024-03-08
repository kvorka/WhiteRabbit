submodule(SphericalHarmonics) vcsum_grid_2
  implicit none; contains
  
  module pure subroutine grid_op_2_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(4, sumNS, grid)
    
    gin(1:2,1:2,1:this%nFourier) => grid(1:4*this%nFourier)
    gout(1:2,1:this%nFourier)    => grid(1:2*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 2
        gout(i2,i) = gin(1,i2,i) * gin(2,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(2, grid, sumNS)
    
  end subroutine grid_op_2_vcsum_sub
  
end submodule vcsum_grid_2