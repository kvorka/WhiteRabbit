submodule (SphericalHarmonics) vcst_grid_2
  implicit none; contains
  
  module pure subroutine grid_op_2_vcst_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl)                        :: tmp
    real(kind=dbl), pointer               :: gin(:,:,:), gout(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(12, sumNS, grid)
    
    gin(1:6,1:2,1:this%nFourier)  => grid(1:12*this%nFourier)
    gout(1:5,1:2,1:this%nFourier) => grid(1:10*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 2
        tmp = gin(1,i2,i)
        
        gout(1,i2,i) = gin(2,i2,i) * tmp
        gout(2,i2,i) = gin(3,i2,i) * tmp
        gout(3,i2,i) = gin(4,i2,i) * tmp
        gout(4,i2,i) = gin(5,i2,i) * tmp
        gout(5,i2,i) = gin(6,i2,i) * tmp
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(10, grid, sumNS)
    
  end subroutine grid_op_2_vcst_sub
  
end submodule vcst_grid_2