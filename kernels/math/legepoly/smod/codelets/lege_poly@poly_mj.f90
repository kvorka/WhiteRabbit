submodule (lege_poly) poly_mj
  implicit none; contains
  
  module procedure mjrec_sub
    integer :: i2
    
    !$omp simd
    do i2 = 1, 16
      pmj2(i2) = cff(3) * pmj1(i2)
      pmj1(i2) = pmj(i2)
      pmj(i2)  = ( cff(1) * cosx2(i2) - cff(2) ) * pmj(i2) - pmj2(i2)
    end do
    
  end procedure mjrec_sub
  
end submodule poly_mj