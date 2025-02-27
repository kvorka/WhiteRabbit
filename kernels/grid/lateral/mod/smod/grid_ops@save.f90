submodule (grid_ops) save
  implicit none; contains
  
  module procedure grid_op_save_sub
    integer :: i2, i1
    
    do i2 = 0, nfour-1
      !$omp simd
      do i1 = 1, step
        grid(i1+i2*padding) = tempgrid(i1,i2+1)
      end do
    end do
    
  end procedure grid_op_save_sub
  
end submodule save