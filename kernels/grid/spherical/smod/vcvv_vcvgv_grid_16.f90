submodule (SphericalHarmonics) vcvv_vcvgv_grid_16
  implicit none; contains
  
  module pure subroutine grid_op_16_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:,:)
    real(kind=dbl), allocatable           :: tmp(:)
    
    call this%fourtrans%exec_c2r_sub( 240, sumNS(1), grid(1) )
    
    allocate( tmp(3) )
    
    gin(1:3,1:5,1:16,1:this%nFourier) => grid(1:240*this%nFourier)
    gout(1:4,1:16,1:this%nFourier)    => grid(1: 64*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, 16
        tmp = gin(1:3,1,i2,i)
        
        gout(1,i2,i) = sum( gin(1:3,2,i2,i) * tmp )
        gout(2,i2,i) = sum( gin(1:3,3,i2,i) * tmp )
        gout(3,i2,i) = sum( gin(1:3,4,i2,i) * tmp )
        gout(4,i2,i) = sum( gin(1:3,5,i2,i) * tmp )
      end do
    end do
    
    deallocate( tmp )
    
    call this%fourtrans%exec_r2c_sub( 64, grid(1), sumNS(1) )
    
  end subroutine grid_op_16_vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv_grid_16