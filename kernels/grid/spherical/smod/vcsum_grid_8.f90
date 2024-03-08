submodule(SphericalHarmonics) vcsum_grid_8
  implicit none; contains
  
  module pure subroutine grid_op_8_vcsum_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(16, sumNS, grid)
    
    gin(1:2,1:8,1:this%nFourier) => grid(1:16*this%nFourier)
    gout(1:8,1:this%nFourier)    => grid(1:8*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 8
        gout(i2,i) = gin(1,i2,i) * gin(2,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(8, grid, sumNS)
    
  end subroutine grid_op_8_vcsum_sub
  
end submodule vcsum_grid_8