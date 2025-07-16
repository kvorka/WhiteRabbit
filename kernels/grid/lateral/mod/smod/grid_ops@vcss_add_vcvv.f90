submodule (grid_ops) vcss_add_vcvv
  implicit none; contains
  
  module procedure grid_op_vcss_add_vcvv_sub
    integer                 :: i1, i2, i3
    real(kind=dbl), pointer :: gout(:,:), gin(:,:,:), gtmp(:,:)
    
    gin(1:16,1:8,1:nfour) => grid(:,1:8*nfour)
    gout(1:16,1:nfour)    => grid(:,1:  nfour)
    gtmp(1:16,1:8)        => tempgrid(:,1:8)
    
    do i3 = 1, nfour
      do i2 = 1, 8
        !$omp simd
        do i1 = 1, 16
          gtmp(i1,i2) = gin(i1,i2,i3)
        end do
      end do
      
      !$omp simd
      do i1 = 1, 16
        gout(i1,i3) = gtmp(i1,7) * gtmp(i1,8)
      end do
      
      do i2 = 1, 3
        !$omp simd
        do i1 = 1, 16
          gout(i1,i3) = gout(i1,i3) + gtmp(i1,i2) * gtmp(i1,i2+3)
        end do
      end do
    end do
    
    gin  => null()
    gout => null()
    gtmp => null()
    
  end procedure grid_op_vcss_add_vcvv_sub
  
end submodule vcss_add_vcvv