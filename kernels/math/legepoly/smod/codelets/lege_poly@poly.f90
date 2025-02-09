submodule (lege_poly) poly
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
  
  module procedure mjrec_sub
    integer :: i2
    
    !$omp simd
    do i2 = 1, step
      pmj2(i2) = cff(3) * pmj1(i2)
      pmj1(i2) = pmj(i2)
      pmj(i2)  = ( cff(1) * cosx2(i2) - cff(2) ) * pmj(i2) - pmj2(i2)
    end do
    
  end procedure mjrec_sub
  
end submodule poly