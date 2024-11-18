submodule (lege_poly) step8p
  implicit none; contains
  
  module procedure mmset_8_sub
    integer :: i2
    
    select case (mj)
      case (1)
        do i2 = 1, 8
          pmm(i2) = this%abmj(1,1)
        end do
      
      case default
        do i2 = 1, 8
          pmm(i2) = this%abmj(1,mj) * sinx(i2) * pmm(i2)
        end do
    end select
    
    do i2 = 1, 8
      pmj(i2,1) = zero
      pmj(i2,2) = zero
      pmj(i2,3) = pmm(i2)
    end do
    
  end procedure mmset_8_sub
  
  module procedure recursion_8_sub
    integer :: i2
    
    do i2 = 1, 8
      pmj(i2,1) = pmj(i2,2)
      pmj(i2,2) = pmj(i2,3)
      pmj(i2,3) = cosx(i2) * this%abmj(2,mj) * pmj(i2,2) - this%abmj(1,mj) * pmj(i2,1)
    end do
    
  end procedure recursion_8_sub
  
  module procedure recursion2_8_sub
    integer :: i2
    
    do i2 = 1, 8
      pmj(i2,1) = pmj(i2,2)
      pmj(i2,2) = pmj(i2,3)
      pmj(i2,3) = cosx(i2) * this%abmj(2,mj) * pmj(i2,2) - this%abmj(1,mj) * pmj(i2,1)
      
      pmj(i2,1) = pmj(i2,2)
      pmj(i2,2) = pmj(i2,3)
      pmj(i2,3) = cosx(i2) * this%abmj(2,mj+1) * pmj(i2,2) - this%abmj(1,mj+1) * pmj(i2,1)
    end do
    
  end procedure recursion2_8_sub
  
end submodule step8p