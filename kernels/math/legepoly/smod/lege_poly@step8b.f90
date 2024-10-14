submodule (lege_poly) step8b
  implicit none; contains
    
  module pure subroutine backward_sum_8_sub(this, nb, legep, cc, legesum)
    class(T_legep),    intent(in)    :: this
    integer,           intent(in)    :: nb
    real(kind=dbl),    intent(in)    :: legep(8)
    complex(kind=dbl), intent(in)    :: cc(nb)
    complex(kind=dbl), intent(inout) :: legesum(8,nb)
    integer                          :: i1, i2
    
    do concurrent ( i1 = 1:nb, i2 = 1:8 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * legep(i2)
    end do
    
  end subroutine backward_sum_8_sub
  
  module pure subroutine backward_rcb_8_sub(this, nb, sumsym, sumasym, sumN, sumS)
    class(T_legep),    intent(in)    :: this
    integer,           intent(in)    :: nb
    complex(kind=dbl), intent(in)    :: sumsym(8,nb), sumasym(8,nb)
    complex(kind=dbl), intent(inout) :: sumN(nb,8), sumS(nb,8)
    
    sumN = transpose( sumsym + sumasym )
    sumS = transpose( sumsym - sumasym )
    
  end subroutine backward_rcb_8_sub
  
  module pure subroutine backward_legesum_8_sub(this, it, nb, cc, sumN, sumS)
    class(T_legep),    intent(in)  :: this
    integer,           intent(in)  :: it, nb
    complex(kind=dbl), intent(in)  :: cc(nb,*)
    complex(kind=dbl), intent(out) :: sumN(*), sumS(*)
    integer                        :: m, j, mj, i2
    real(kind=dbl),    allocatable :: pmj2(:), pmj1(:), pmj0(:), pmm(:), csx(:), snx(:)
    complex(kind=dbl), allocatable :: ssm(:), asm(:)
    
    allocate( pmj2(8), pmj1(8), pmj0(8), pmm(8), csx(8), snx(8), ssm(8*nb), asm(8*nb) )
    
    do concurrent ( i2 = 0:7 )
      csx(i2+1) = this%roots(it+i2)
      snx(i2+1) = sqrt(1-this%roots(it+i2)**2)
    end do
    
    call zero_carray_sub( 8*nb*(this%jmax+1), sumN(1) )
    call zero_carray_sub( 8*nb*(this%jmax+1), sumS(1) )
    
    do m = 0, this%jmax
      call zero_carray_sub( 8*nb, ssm(1) )
      call zero_carray_sub( 8*nb, asm(1) )
      
      !j = m
        mj = m*(this%jmax+1)-(m-2)*(m+1)/2
        
        call this%mmset_8_sub( m, snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_8_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
        
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        
        call this%recursion_8_sub( mj-1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_8_sub( nb, pmj0(1), cc(1,mj-1), asm(1) )
        
        call this%recursion_8_sub( mj, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_8_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        call this%recursion_8_sub( mj+1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_8_sub( nb, pmj0(1), cc(1,mj+1), asm(1) )
      end if
      
      call this%backward_rcb_8_sub( nb, ssm(1), asm(1), sumN(8*nb*m+1), sumS(8*nb*m+1) )
    end do
    
    deallocate( pmj2, pmj1, pmj0, pmm, csx, snx, asm, ssm )
    
  end subroutine backward_legesum_8_sub
  
end submodule step8b