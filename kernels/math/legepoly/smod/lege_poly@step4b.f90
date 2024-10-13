submodule (lege_poly) step4b
  implicit none; contains
  
  module pure subroutine backward_sum_4_sub(this, nb, legep, cc, legesum)
    class(T_legep),    intent(in)    :: this
    integer,           intent(in)    :: nb
    real(kind=dbl),    intent(in)    :: legep(4)
    complex(kind=dbl), intent(in)    :: cc(nb)
    complex(kind=dbl), intent(inout) :: legesum(4,nb)
    integer                          :: i1, i2
    
    do concurrent ( i1 = 1:nb, i2 = 1:4 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * legep(i2)
    end do
    
  end subroutine backward_sum_4_sub
  
  module pure subroutine backward_rcb_4_sub(this, nb, sumsym, sumasym, sumN, sumS)
    class(T_legep),    intent(in)    :: this
    integer,           intent(in)    :: nb
    complex(kind=dbl), intent(in)    :: sumsym(4,nb), sumasym(4,nb)
    complex(kind=dbl), intent(inout) :: sumN(nb,4), sumS(nb,4)
    
    sumN = transpose( sumsym + sumasym )
    sumS = transpose( sumsym - sumasym )
    
  end subroutine backward_rcb_4_sub
  
  module pure subroutine backward_legesum_4_sub(this, nb, roots, cc, sumN, sumS)
    class(T_legep),    intent(in)  :: this
    integer,           intent(in)  :: nb
    real(kind=dbl),    intent(in)  :: roots(4)
    complex(kind=dbl), intent(in)  :: cc(nb,*)
    complex(kind=dbl), intent(out) :: sumN(*), sumS(*)
    integer                        :: m, j, mj, i2
    real(kind=dbl),    allocatable :: pmj2(:), pmj1(:), pmj0(:), pmm(:), csx(:), snx(:)
    complex(kind=dbl), allocatable :: ssm(:), asm(:)
    
    allocate( pmj2(4), pmj1(4), pmj0(4), pmm(4), csx(4), snx(4), ssm(4*nb), asm(4*nb) )
    
    do concurrent ( i2 = 1:4 )
      csx(i2) = roots(i2)
      snx(i2) = sqrt(1-csx(i2)**2)
    end do
    
    call zero_carray_sub( 4*nb*(this%jmax+1), sumN(1) )
    call zero_carray_sub( 4*nb*(this%jmax+1), sumS(1) )
    
    do m = 0, this%jmax
      call zero_carray_sub( 4*nb, ssm(1) )
      call zero_carray_sub( 4*nb, asm(1) )
      
      !j = m
        mj = m*(this%jmax+1)-(m-2)*(m+1)/2
        
        call this%mmset_4_sub( m, snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_4_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
        
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        
        call this%recursion_4_sub( mj-1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_4_sub( nb, pmj0(1), cc(1,mj-1), asm(1) )
        
        call this%recursion_4_sub( mj, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_4_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        call this%recursion_4_sub( mj+1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%backward_sum_4_sub( nb, pmj0(1), cc(1,mj+1), asm(1) )
      end if
      
      call this%backward_rcb_4_sub( nb, ssm(1), asm(1), sumN(4*nb*m+1), sumS(4*nb*m+1) )
    end do
    
    deallocate( pmj2, pmj1, pmj0, pmm, csx, snx, asm, ssm )
    
  end subroutine backward_legesum_4_sub
  
end submodule step4b