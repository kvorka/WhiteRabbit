submodule (grid_ops) vcss
  implicit none; contains
  
  module procedure grid_op_vcss_sub
    integer                 :: i1, i2, i3
    real(kind=dbl), pointer :: gout(:,:), gin(:,:,:), gtmp(:,:)
    
    gin(1:16,1:2,1:nfour) => grid(:,1:2*nfour)
    gout(1:16,1:nfour)    => grid(:,1:  nfour)
    gtmp(1:16,1:2)        => tempgrid(:,1:2)
    
    do i3 = 1, nfour
      do i2 = 1, 2
        !$omp simd
        do i1 = 1, 16
          gtmp(i1,i2) = gin(i1,i2,i3)
        end do
      end do
      
      !$omp simd
      do i1 = 1, 16
        gout(i1,i3) = gtmp(i1,1) * gtmp(i1,2)
      end do
    end do
    
    gin  => null()
    gout => null()
    gtmp => null()
    
  end procedure grid_op_vcss_sub
  
end submodule vcss