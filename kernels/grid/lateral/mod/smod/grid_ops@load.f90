submodule (grid_ops) load
  implicit none; contains
  
  module procedure grid_op_load_sub
    integer :: i2, i1
    
    do i2 = 0, nfour-1
      !$omp simd
      do i1 = 1, step
        tempgrid(i1,i2+1) = grid(i1+i2*padding)
      end do
    end do
    
  end procedure grid_op_load_sub
  
end submodule load