submodule (grid_ops) vcvv
  implicit none; contains
  
  module procedure grid_op_vcvv_sub
    integer                     :: i1, i2
    real(kind=dbl), pointer     :: gout(:,:), gin(:,:,:,:)
    real(kind=dbl), allocatable :: tmp1(:), tmp2(:)
    
    gin(1:step,1:3,1:2,1:nfour) => grid(:,1:6*nfour)
    gout(1:step,1:nfour)        => grid(:,1:  nfour)
    
    allocate( tmp1(step), tmp2(step) )
    
    do i1 = 1, nfour
      tmp1 = gin(1:step,1,1,i1)
      tmp2 = gin(1:step,1,2,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,i1) = tmp1(i2) * tmp2(i2)
      end do
      
      tmp1 = gin(1:step,2,1,i1)
      tmp2 = gin(1:step,2,2,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,i1) = gout(i2,i1) + tmp1(i2) * tmp2(i2)
      end do
      
      tmp1 = gin(1:step,3,1,i1)
      tmp2 = gin(1:step,3,2,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,i1) = gout(i2,i1) + tmp1(i2) * tmp2(i2)
      end do
    end do
    
    deallocate( tmp1, tmp2 )
    
    gin  => null()
    gout => null()
    
  end procedure grid_op_vcvv_sub
  
end submodule vcvv