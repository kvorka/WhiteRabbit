submodule (math) nulify
  implicit none; contains
  
  module procedure zero_rarray_sub
    integer :: i
    
    !$omp simd
    do i = 1, length
      arr(i) = zero
    end do
    
  end procedure zero_rarray_sub
  
  module procedure zero_carray_sub
    integer :: i
    
    !$omp simd
    do i = 1, length
      arr(i) = czero
    end do
    
  end procedure zero_carray_sub
  
end submodule nulify