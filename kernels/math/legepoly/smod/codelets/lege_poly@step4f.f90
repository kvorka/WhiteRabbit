submodule (lege_poly) step4f
  implicit none; contains
  
  module procedure forward_sum_4_sub
    integer :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + sum( legep(1:4) * legesum(1:4,i1) )
    end do
    
  end procedure forward_sum_4_sub
  
  module procedure forward_rcb_4_sub
    integer :: i1, i2
    
    do concurrent ( i2 = 1:nf, i1 = 1:4 )
      sumsym(i1,i2)  = cmplx( sumN(i1,i2,1) + sumS(i1,i2,1), sumN(i1,i2,2) + sumS(i1,i2,2), kind=dbl )
      sumasym(i1,i2) = cmplx( sumN(i1,i2,1) - sumS(i1,i2,1), sumN(i1,i2,2) - sumS(i1,i2,2), kind=dbl )
    end do
    
    do concurrent ( i2 = 1:nf, i1 = 1:4 )
      sumsym(i1,i2)  = w(i1) * sumsym(i1,i2)
      sumasym(i1,i2) = w(i1) * sumasym(i1,i2)
    end do
    
  end procedure forward_rcb_4_sub
  
  module procedure forward_legesum_4_sub
    integer                        :: j, m, mj, i2
    real(kind=dbl),    allocatable :: pmj2(:), pmj1(:), pmj0(:), pmm(:)
    complex(kind=dbl), allocatable :: ssm(:), asm(:)
    
    allocate( pmj2(4), pmj1(4), pmj0(4), pmm(4), ssm(4*nf), asm(4*nf) )
    
    do m = 0, this%jmax
      call this%forward_rcb_4_sub( nf, rwork(1,3), sumN(1,1,1,m), sumS(1,1,1,m), ssm(1), asm(1) )
      
      !j = m
        mj = m*(this%jmax+1)-(m-2)*(m+1)/2
        
        call this%mmset_4_sub( mj, rwork(1,2), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_4_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
      
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        
        call this%recursion_4_sub( mj-1, rwork(1,1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_4_sub( nf, pmj0(1), asm(1), cr(1,mj-1) )
        
        call this%recursion_4_sub( mj, rwork(1,1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_4_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
      end do
      
      if ( mod(this%jmax-m,2) /= 0 ) then
        call this%recursion_4_sub( mj+1, rwork(1,1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_4_sub( nf, pmj0(1), asm(1), cr(1,mj+1) )
      end if
    end do
    
    deallocate( pmj2, pmj1, pmj0, pmm, asm, ssm )
    
  end procedure forward_legesum_4_sub
  
end submodule step4f