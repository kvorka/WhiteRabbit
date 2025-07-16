submodule (grid_ops) vcvxv
  use math
  implicit none; contains
  
  module procedure grid_op_vcvxv_sub
    integer                     :: i1, i2, i3, i4
    real(kind=dbl), pointer     :: gout(:,:,:), gin(:,:,:,:), gtmp(:,:,:)
    
    gin(1:16,1:3,1:2,1:nfour) => grid(:,1:6*nfour)
    gout(1:16,1:3,1:nfour)    => grid(:,1:3*nfour)
    gtmp(1:16,1:3,1:2)        => tempgrid(:,1:6)
    
    do i4 = 1, nfour
      do i3 = 1, 2
        do i2 = 1, 3
          !$omp simd
          do i1 = 1, 16
            gtmp(i1,i2,i3) = gin(i1,i2,i3,i4)
          end do
        end do
      end do
      
      !$omp simd
      do i1 = 1, 16
        gout(i1,1,i4) = gtmp(i1,2,1) * gtmp(i1,3,2) - gtmp(i1,3,1) * gtmp(i1,2,2)
        gout(i1,2,i4) = gtmp(i1,3,1) * gtmp(i1,1,2) - gtmp(i1,1,1) * gtmp(i1,3,2)
        gout(i1,3,i4) = gtmp(i1,1,1) * gtmp(i1,2,2) - gtmp(i1,2,1) * gtmp(i1,1,2)
      end do
    end do
    
    gin  => null()
    gout => null()
    gtmp => null()
    
  end procedure grid_op_vcvxv_sub
  
end submodule vcvxv