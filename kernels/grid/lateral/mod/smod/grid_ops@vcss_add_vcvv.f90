submodule (grid_ops) vcss_add_vcvv
  implicit none; contains
  
  module pure subroutine grid_op_vcss_add_vcvv_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(step,*)
    integer                               :: i1, i2
    real(kind=dbl), allocatable           :: tmp(:), tmp1(:)
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    gin(1:step,1:8,1:nfour) => grid(:,1:8*nfour)
    gout(1:step,1:nfour)    => grid(:,1:  nfour)
    
    allocate( tmp(step), tmp1(step) )
    
    do i1 = 1, nfour
      tmp  = gin(1:step,1,i1)
      tmp1 = gin(1:step,4,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,i1) = tmp(i2) * tmp1(i2)
      end do
      
      tmp  = gin(1:step,2,i1)
      tmp1 = gin(1:step,5,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,i1) = gout(i2,i1) + tmp(i2) * tmp1(i2)
      end do
      
      tmp  = gin(1:step,3,i1)
      tmp1 = gin(1:step,6,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,i1) = gout(i2,i1) + tmp(i2) * tmp1(i2)
      end do
      
      tmp  = gin(1:step,7,i1)
      tmp1 = gin(1:step,8,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,i1) = gout(i2,i1) + tmp(i2) * tmp1(i2)
      end do
    end do
    
    deallocate( tmp, tmp1 )
    
  end subroutine grid_op_vcss_add_vcvv_sub
  
end submodule vcss_add_vcvv