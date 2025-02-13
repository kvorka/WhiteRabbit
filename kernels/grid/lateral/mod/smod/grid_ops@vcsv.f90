submodule (grid_ops) vcsv
  implicit none; contains
  
  module procedure grid_op_vcsv_sub
    integer                 :: i1, i2, i3
    real(kind=dbl), pointer :: gout(:,:,:), gin(:,:,:), gtmp(:,:)
    
    gin(1:step,1:4,1:nfour)  => grid(:,1:4*nfour)
    gout(1:step,1:3,1:nfour) => grid(:,1:3*nfour)
    gtmp(1:step,1:4)         => grid(:,1:4)
    
    do i3 = 1, nfour
      !$omp simd collapse (2)
      do i2 = 1, 4
        do i1 = 1, step
          gtmp(i1,i2) = gin(i1,i2,i3)
        end do
      end do
      
      !$omp simd collapse (2)
      do i2 = 1, 3
        do i1 = 1, step
          gout(i1,i2,i3) = gtmp(i1,1) * gtmp(i1,i2+1)
        end do
      end do
    end do
    
    gin  => null()
    gout => null()
    gtmp => null()
    
  end procedure grid_op_vcsv_sub
  
end submodule vcsv