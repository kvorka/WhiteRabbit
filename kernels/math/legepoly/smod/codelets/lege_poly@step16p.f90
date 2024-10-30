submodule (lege_poly) step16p
  implicit none; contains
  
  module procedure mmset_16_sub
    integer :: i2
    
    if ( m /= 0 ) then
      do concurrent ( i2 = 1:16 )
        pmm(i2) = this%cmm(m) * sinx(i2) * pmm(i2)
      end do
    else
      do concurrent ( i2 = 1:16 )
        pmm(i2) = this%cmm(0)
      end do
    end if
    
    do concurrent ( i2 = 1:16 )
      pmj2(i2) = zero
      pmj1(i2) = zero
      pmj0(i2) = pmm(i2)
    end do
    
  end procedure mmset_16_sub
  
  module procedure recursion_16_sub
    integer :: i2
    
    do concurrent ( i2=1:16 )
      pmj2(i2) = this%bmj(mj) * pmj1(i2)
      pmj1(i2) = this%amj(mj) * pmj0(i2)
      pmj0(i2) = cosx(i2) * pmj1(i2) - pmj2(i2)
    end do
    
  end procedure recursion_16_sub
  
end submodule step16p