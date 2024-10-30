submodule (grid_ops) vcst
  implicit none; contains
  
  module procedure grid_op_vcst_sub
    integer                     :: i1, i2, i3
    real(kind=dbl), pointer     :: gin(:,:,:), gout(:,:,:)
    real(kind=dbl), allocatable :: tmp(:), tmp1(:)
    
    gin(1:step,1:6,1:nfour)  => grid(:,1:6*nfour)
    gout(1:step,1:5,1:nfour) => grid(:,1:5*nfour)
    
    allocate( tmp(step), tmp1(step) )
    
    do i1 = 1, nfour
      tmp = gin(1:step,1,i1)
      
      do i3 = 1, 5
        tmp1 = gin(1:step,i3+1,i1)
        
        do concurrent ( i2 = 1:step )
          gout(i2,i3,i1) = tmp(i2) * tmp1(i2)
        end do
      end do
    end do
    
    deallocate( tmp, tmp1 )
    
    gin  => null()
    gout => null()
    
  end procedure grid_op_vcst_sub
  
end submodule vcst