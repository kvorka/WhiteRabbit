submodule (SphericalHarmonics) vcvv_vcvgv_grid_4
  implicit none; contains
  
  module pure subroutine grid_op_4_vcvv_vcvgv_sub(this, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i
    real(kind=dbl), pointer               :: gin(:), gout(:)
    real(kind=dbl), allocatable           :: tmp(:)
    
    call this%fourtrans%exec_c2r_sub(60, sumNS, grid)
    
    allocate( tmp(3) )
    
    do i = 0, this%nFourier-1
      gin  => grid(1+60*i:60+60*i)
      gout => grid(1+16*i:16+16*i)
      
      tmp = gin(1:3)
        gout(1) = sum( gin( 4: 6) * tmp )
        gout(2) = sum( gin( 7: 9) * tmp )
        gout(3) = sum( gin(10:12) * tmp )
        gout(4) = sum( gin(13:15) * tmp )
      
      tmp = gin(16:18)
        gout(5) = sum( gin(19:21) * tmp )
        gout(6) = sum( gin(22:24) * tmp )
        gout(7) = sum( gin(25:27) * tmp )
        gout(8) = sum( gin(28:30) * tmp )
      
      tmp = gin(31:33)
        gout( 9) = sum( gin(34:36) * tmp )
        gout(10) = sum( gin(37:39) * tmp )
        gout(11) = sum( gin(40:42) * tmp )
        gout(12) = sum( gin(43:45) * tmp )
      
      tmp = gin(46:48)
        gout(13) = sum( gin(49:51) * tmp )
        gout(14) = sum( gin(52:54) * tmp )
        gout(15) = sum( gin(55:57) * tmp )
        gout(16) = sum( gin(58:60) * tmp )
    end do
    
    deallocate( tmp )
    
    call this%fourtrans%exec_r2c_sub(16, grid, sumNS)
    
  end subroutine grid_op_4_vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv_grid_4