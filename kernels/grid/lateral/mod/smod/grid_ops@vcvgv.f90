submodule (grid_ops) vcvgv
  implicit none; contains
  
  module procedure grid_op_vcvgv_sub
    integer                 :: i1, i2, i3, i4
    real(kind=dbl), pointer :: gout(:,:,:), gin(:,:,:,:), gtmp(:,:,:)
    
    gin(1:16,1:3,1:4,1:nfour) => grid(:,1:12*nfour)
    gout(1:16,1:3,1:nfour)    => grid(:,1: 3*nfour)
    gtmp(1:16,1:3,1:4)        => tempgrid(:,1:12)
    
    do i4 = 1, nfour
      do i3 = 1, 4
        do i2 = 1, 3
          !$omp simd
          do i1 = 1, 16
            gtmp(i1,i2,i3) = gin(i1,i2,i3,i4)
          end do
        end do
      end do
      
      do i3 = 1, 3
        !$omp simd
        do i1 = 1, 16
          gout(i1,i3,i4) = zero
        end do
      end do
      
      do i3 = 1, 3
        do i2 = 1, 3
          !$omp simd
          do i1 = 1, 16
            gout(i1,i3,i4) = gout(i1,i3,i4) + gtmp(i1,i2,1) * gtmp(i1,i2,i3+1)
          end do
        end do
      end do
    end do
    
    gin  => null()
    gout => null()
    gtmp => null()
    
  end procedure grid_op_vcvgv_sub
  
end submodule vcvgv