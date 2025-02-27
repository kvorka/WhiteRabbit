submodule (lege_poly) is_rescale
  implicit none; contains
  
  module procedure is_rescale_sub
    integer :: i1, i2
    
    do i2 = 1, this%nrma
      !$omp simd
      do i1 = 1, 4*ncab
        rcab(i1,i2) = this%amj(i2) * rcab(i1,i2)
      end do
    end do
    
  end procedure is_rescale_sub
  
end submodule is_rescale