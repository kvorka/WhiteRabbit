submodule (lege_poly) step16
  implicit none; contains
  
  module procedure backward_sum_16_sub
    integer :: i1, i2
    
    do concurrent ( i1 = 1:nb, i2 = 1:16 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * legep(i2)
    end do
    
  end procedure backward_sum_16_sub
  
  module procedure backward_rcb_16_sub
    integer :: i1, i2
    
    do concurrent ( i2 = 1:nb, i1 = 1:16 )
      sumN(i1,i2,1) = sumsym(i1,i2)%re + sumasym(i1,i2)%re
      sumN(i1,i2,2) = sumsym(i1,i2)%im + sumasym(i1,i2)%im
      sumS(i1,i2,1) = sumsym(i1,i2)%re - sumasym(i1,i2)%re
      sumS(i1,i2,2) = sumsym(i1,i2)%im - sumasym(i1,i2)%im
    end do
    
  end procedure backward_rcb_16_sub
  
  module procedure backward_legesum_16_sub
    integer                        :: m, j, mj, i2
    real(kind=dbl),    allocatable :: pmj2(:), pmj1(:), pmj0(:), pmm(:), csx(:), snx(:)
    complex(kind=dbl), allocatable :: ssm(:), asm(:)
    
    allocate( pmj2(16), pmj1(16), pmj0(16), pmm(16), csx(16), snx(16), ssm(16*nb), asm(16*nb) )
    
    do concurrent ( i2 = 0:15 )
      csx(i2+1) = this%rootsweights(1,it+i2)
      snx(i2+1) = this%rootsweights(2,it+i2)
    end do
    
    do m = 0, this%jmax
      call zero_carray_sub( 16*nb, ssm(1) )
      call zero_carray_sub( 16*nb, asm(1) )
      
      !j = m
        mj = m*(this%jmax+1)-(m-2)*(m+1)/2
        
        call this%mmset_16_sub( m, snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_16_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
        
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        
        call this%recursion_16_sub( mj-1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_16_sub( nb, pmj0(1), cc(1,mj-1), asm(1) )
        
        call this%recursion_16_sub( mj, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_16_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        call this%recursion_16_sub( mj+1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_16_sub( nb, pmj0(1), cc(1,mj+1), asm(1) )
      end if
      
      call this%backward_rcb_16_sub( nb, ssm(1), asm(1), sumN(1,1,1,m), sumS(1,1,1,m) )
    end do
    
    deallocate( pmj2, pmj1, pmj0, pmm, csx, snx, asm, ssm )
    
  end procedure backward_legesum_16_sub

end submodule step16