submodule (grid_ops) vcvxv
  use math
  implicit none; contains
  
  module pure subroutine grid_op_vcvxv_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(step,*)
    integer                               :: i1, i2
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:,:)
    real(kind=dbl), allocatable           :: tmp11(:), tmp12(:), tmp21(:), tmp22(:)
    
    gin(1:step,1:3,1:2,1:nfour) => grid(:,1:6*nfour)
    gout(1:step,1:3,1:nfour)    => grid(:,1:3*nfour)
    
    allocate( tmp11(step), tmp12(step), tmp21(step), tmp22(step) )
    
    do i1 = 1, nfour
      tmp11 = gin(1:step,2,1,i1)
      tmp12 = gin(1:step,3,1,i1)
      tmp21 = gin(1:step,2,2,i1)
      tmp22 = gin(1:step,3,2,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,1,i1) = tmp11(i2) * tmp22(i2) - tmp12(i2) * tmp21(i2)
      end do
      
      tmp11 = gin(1:step,3,1,i1)
      tmp12 = gin(1:step,1,1,i1)
      tmp21 = gin(1:step,3,2,i1)
      tmp22 = gin(1:step,1,2,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,2,i1) = tmp11(i2) * tmp22(i2) - tmp12(i2) * tmp21(i2)
      end do
      
      tmp11 = gin(1:step,1,1,i1)
      tmp12 = gin(1:step,2,1,i1)
      tmp21 = gin(1:step,1,2,i1)
      tmp22 = gin(1:step,2,2,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,3,i1) = tmp11(i2) * tmp22(i2) - tmp12(i2) * tmp21(i2)
      end do
    end do
    
    deallocate( tmp11, tmp12, tmp21, tmp22 )
    
    gin  => null()
    gout => null()
    
  end subroutine grid_op_vcvxv_sub
  
end submodule vcvxv