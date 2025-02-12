submodule (grid_ops) vcst
  implicit none; contains
  
  module procedure grid_op_vcst_sub
    integer                 :: i1, i2, i3
    real(kind=dbl), pointer :: gin(:,:,:), gout(:,:,:), gtmp(:,:)
    
    gin(1:step,1:6,1:nfour)  => grid(:,1:6*nfour)
    gout(1:step,1:5,1:nfour) => grid(:,1:5*nfour)
    gtmp(1:step,1:6)         => tempgrid(:,1:6)
    
    do i3 = 1, nfour
      !$omp simd collapse (2)
      do i2 = 1, 6
        do i1 = 1, step
          gtmp(i1,i2) = gin(i1,i2,i3)
        end do
      end do
      
      !$omp simd collapse (2)
      do i2 = 1, 5
        do i1 = 1, step
          gout(i1,i2,i3) = gtmp(i1,1) * gtmp(i1,i2+1)
        end do
      end do
    end do
    
    gin  => null()
    gout => null()
    gtmp => null()
    
  end procedure grid_op_vcst_sub
  
end submodule vcst