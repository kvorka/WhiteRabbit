submodule (lege_poly) poly_mm
  implicit none; contains
  
  module procedure mmset_sub
    integer :: i2
    
    select case (ma)
      case (1)
        !$omp simd
        do i2 = 1, step
          pmm(i2) = cff
        end do
      
      case default
        !$omp simd
        do i2 = 1, step
          pmm(i2) = cff * sinx(i2) * pmm(i2)
        end do
    end select
    
    !$omp simd
    do i2 = 1, step
      pmj2(i2) = zero
      pmj1(i2) = zero
      pmj(i2)  = pmm(i2) / cosx(i2)
    end do
    
  end procedure mmset_sub
  
end submodule poly_mm