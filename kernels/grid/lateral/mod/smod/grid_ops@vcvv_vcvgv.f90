submodule (grid_ops) vcvv_vcvgv
  implicit none; contains
  
  module pure subroutine grid_op_vcvv_vcvgv_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(step,*)
    integer                               :: i1, i2, i3
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:,:)
    real(kind=dbl), allocatable           :: tmp11(:), tmp12(:), tmp13(:), tmp(:)
    
    allocate( tmp11(step), tmp12(step), tmp13(step), tmp(step) )
    
    gin(1:step,1:3,1:5,1:nfour) => grid(:,1:15*nfour)
    gout(1:step,1:4,1:nfour)    => grid(:,1: 4*nfour)
    
    do i1 = 1, nfour
      tmp11 = gin(1:step,1,1,i1)
      tmp12 = gin(1:step,2,1,i1)
      tmp13 = gin(1:step,3,1,i1)
      
      do i3 = 1, 4
        tmp = gin(1:step,1,i3+1,i1)
          do concurrent ( i2 = 1:step )
            gout(i2,i3,i1) = tmp11(i2) * tmp(i2)
          end do
        
        tmp = gin(1:step,2,i3+1,i1)
          do concurrent ( i2 = 1:step )
            gout(i2,i3,i1) = gout(i2,i3,i1) + tmp12(i2) * tmp(i2)
          end do
        
        tmp = gin(1:step,3,i3+1,i1)
          do concurrent ( i2 = 1:step )
            gout(i2,i3,i1) = gout(i2,i3,i1) + tmp13(i2) * tmp(i2)
          end do
      end do
    end do
    
    deallocate( tmp11, tmp12, tmp13, tmp )
    
  end subroutine grid_op_vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv