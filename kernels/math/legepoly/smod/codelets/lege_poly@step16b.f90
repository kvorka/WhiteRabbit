submodule (lege_poly) step16
  implicit none; contains
  
  module procedure backward_sum_16_sub
    integer :: i1, i2
    
    do i1 = 1, nb
      do i2 = 1, 16
        legesum(i2,i1) = legesum(i2,i1) + cc(i1) * legep(i2)
      end do
    end do
    
  end procedure backward_sum_16_sub
  
  module procedure backward_sum2_16_sub
    integer :: i1, i2, i3
    
    do i3 = 1, 2
      do i2 = 1, nb
        do i1 = 1, 16
          legesum(i1,i2,i3) = legesum(i1,i2,i3) + legep(i1,i3) * cc(i2,i3)
        end do
      end do
    end do
    
  end procedure backward_sum2_16_sub
  
  module procedure backward_rcb_16_sub
    integer :: i1, i2
    
    do i2 = 1, nb
      do i1 = 1, 16
        sumN(i1,i2,1) = swork(i1,i2,2)%re + swork(i1,i2,1)%re
        sumS(i1,i2,1) = swork(i1,i2,2)%re - swork(i1,i2,1)%re
        sumN(i1,i2,2) = swork(i1,i2,2)%im + swork(i1,i2,1)%im
        sumS(i1,i2,2) = swork(i1,i2,2)%im - swork(i1,i2,1)%im
      end do
    end do
    
  end procedure backward_rcb_16_sub
  
  module procedure backward_legesum_16_sub
    integer                        :: m, j, mj
    real(kind=dbl),    allocatable :: pmm(:), pmj(:)
    complex(kind=dbl), allocatable :: swork(:)
    
    allocate( pmm(16), pmj(48), swork(32*nb) )
    
    do m = 0, this%jmax
      call zero_carray_sub( 32*nb, swork(1) )
      
      !j = m
        mj = m*(this%jmax+1)-(m-2)*(m+1)/2
        
        call this%mmset_16_sub( mj, rw(1,2), pmm(1), pmj(1) )
        call this%backward_sum_16_sub( nb, pmj(33), cc(1,mj), swork(16*nb+1) )
        
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        
        call this%recursion2_16_sub( mj-1, rw(1,1), pmj(1) )
        call this%backward_sum2_16_sub( nb, pmj(17), cc(1,mj-1), swork(1) )
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        call this%recursion_16_sub( mj+1, rw(1,1), pmj(1) )
        call this%backward_sum_16_sub( nb, pmj(33), cc(1,mj+1), swork(1) )
      end if
      
      call this%backward_rcb_16_sub( nb, swork(1), sumN(1,m), sumS(1,m) )
    end do
    
    deallocate( pmm, pmj, swork )
    
  end procedure backward_legesum_16_sub

end submodule step16