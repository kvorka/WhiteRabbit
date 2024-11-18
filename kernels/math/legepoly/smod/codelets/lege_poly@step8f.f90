submodule (lege_poly) step8f
  implicit none; contains
  
  module procedure forward_sum_8_sub
    integer :: i1
    
    do concurrent ( i1 = 1:nf )
      cr(i1) = cr(i1) + sum( legep(1:8) * legesum(1:8,i1) )
    end do
    
  end procedure forward_sum_8_sub
  
  module procedure forward_sum2_8_sub
    integer :: i1, i2
    
    do concurrent ( i2 = 1:2, i1 = 1:nf )
      cr(i1,i2) = cr(i1,i2) + sum( legep(1:8,i2) * legesum(1:8,i1,i2) )
    end do
    
  end procedure forward_sum2_8_sub
  
  module procedure forward_rcb_8_sub
    integer :: i1, i2, i3
    
    do concurrent ( i2 = 1:nf, i1 = 1:8 )
      swork(i1,i2,1) = cmplx( sumN(i1,i2,1) - sumS(i1,i2,1), sumN(i1,i2,2) - sumS(i1,i2,2), kind=dbl )
      swork(i1,i2,2) = cmplx( sumN(i1,i2,1) + sumS(i1,i2,1), sumN(i1,i2,2) + sumS(i1,i2,2), kind=dbl )
    end do
    
    do concurrent ( i3 = 1:2, i2 = 1:nf, i1 = 1:8 )
      swork(i1,i2,i3)  = w(i1) * swork(i1,i2,i3)
    end do
    
  end procedure forward_rcb_8_sub
  
  module procedure forward_legesum_8_sub
    integer                        :: j, m, mj
    real(kind=dbl),    allocatable :: pmm(:), pmj(:)
    complex(kind=dbl), allocatable :: swork(:)
    
    allocate( pmm(8), pmj(24), swork(16*nf) )
    
    do m = 0, this%jmax
      call this%forward_rcb_8_sub( nf, rw(1,3), sumN(1,m), sumS(1,m), swork(1) )
      
      !j = m
        mj = m*(this%jmax+1)-(m-2)*(m+1)/2
        
        call this%mmset_8_sub( mj, rw(1,2), pmm(1), pmj(1) )
        call this%forward_sum_8_sub( nf, pmj(17), swork(8*nf+1), cr(1,mj) )
      
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        
        call this%recursion2_8_sub( mj-1, rw(1,1), pmj(1) )
        call this%forward_sum2_8_sub( nf, pmj(9), swork(1), cr(1,mj-1) )
      end do
      
      if ( mod(this%jmax-m,2) /= 0 ) then
        call this%recursion_8_sub( mj+1, rw(1,1), pmj(1) )
        call this%forward_sum_8_sub( nf, pmj(17), swork(1), cr(1,mj+1) )
      end if
    end do
    
    deallocate( pmm, pmj, swork )
    
  end procedure forward_legesum_8_sub
  
end submodule step8f