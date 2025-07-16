submodule (grid_ops) vcsv
  implicit none; contains
  
  module procedure grid_op_vcsv_sub
    integer                 :: i1, i2, i3
    real(kind=dbl), pointer :: gout(:,:,:), gin(:,:,:), gtmp(:,:)
    
    gin(1:16,1:4,1:nfour)  => grid(:,1:4*nfour)
    gout(1:16,1:3,1:nfour) => grid(:,1:3*nfour)
    gtmp(1:16,1:4)         => grid(:,1:4)
    
    do i3 = 1, nfour
      do i2 = 1, 4
        !$omp simd
        do i1 = 1, 16
          gtmp(i1,i2) = gin(i1,i2,i3)
        end do
      end do
      
      do i2 = 1, 3
        !$omp simd
        do i1 = 1, 16
          gout(i1,i2,i3) = gtmp(i1,1) * gtmp(i1,i2+1)
        end do
      end do
    end do
    
    gin  => null()
    gout => null()
    gtmp => null()
    
  end procedure grid_op_vcsv_sub
  
end submodule vcsv