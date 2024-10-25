submodule (grid_ops) vcsv
  implicit none; contains
  
  module pure subroutine grid_op_vcsv_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(step,*)
    integer                               :: i1, i2, i3
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:)
    real(kind=dbl), allocatable           :: tmp(:), tmp1(:)
    
    gin(1:step,1:4,1:nfour)  => grid(:,1:4*nfour)
    gout(1:step,1:3,1:nfour) => grid(:,1:3*nfour)
    
    allocate( tmp(step), tmp1(step) )
    
    do i1 = 1, nfour
      tmp = gin(1:step,1,i1)
      
      do i3 = 1, 3
        tmp1 = gin(1:step,i3+1,i1)
        
        do concurrent ( i2 = 1:step )
          gout(i2,i3,i1) = tmp(i2) * tmp1(i2)
        end do
      end do
    end do
    
    deallocate( tmp, tmp1 )
    
    gin  => null()
    gout => null()
    
  end subroutine grid_op_vcsv_sub
  
end submodule vcsv