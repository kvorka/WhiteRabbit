submodule (lege_poly) step16f
  implicit none; contains
  
  module pure subroutine forward_sum_16_sub(this, nf, legep, legesum, cr)
    class(T_legep),    intent(in)    :: this
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: legep(16)
    complex(kind=dbl), intent(in)    :: legesum(16,nf)
    complex(kind=dbl), intent(inout) :: cr(nf)
    integer                          :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + sum( legep(1:16) * legesum(1:16,i1) )
    end do
    
  end subroutine forward_sum_16_sub
  
  module pure subroutine forward_rcb_16_sub(this, nf, w, sumN, sumS, sumsym, sumasym)
    class(T_legep),    intent(in)  :: this
    integer,           intent(in)  :: nf
    real(kind=dbl),    intent(in)  :: w(16)
    complex(kind=dbl), intent(in)  :: sumN(nf,16), sumS(nf,16)
    complex(kind=dbl), intent(out) :: sumsym(16,nf), sumasym(16,nf)
    integer                        :: i1, i2
    
    sumsym  = transpose( sumN + sumS )
    sumasym = transpose( sumN - sumS )
    
    do concurrent ( i1 = 1:nf, i2 = 1:16 )
      sumsym(i2,i1)  = w(i2) * sumsym(i2,i1)
      sumasym(i2,i1) = w(i2) * sumasym(i2,i1)
    end do
    
  end subroutine forward_rcb_16_sub
  
  module pure subroutine forward_legesum_16_sub(this, nf, roots, weights, sumN, sumS, cr)
    class(T_legep),    intent(in)    :: this
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: roots(16), weights(16)
    complex(kind=dbl), intent(in)    :: sumN(*), sumS(*)
    complex(kind=dbl), intent(inout) :: cr(nf,*)
    integer                          :: j, m, mj, i2
    real(kind=dbl),    allocatable   :: pmj2(:), pmj1(:), pmj0(:), pmm(:), csx(:), snx(:)
    complex(kind=dbl), allocatable   :: ssm(:), asm(:)
    
    allocate( pmj2(16), pmj1(16), pmj0(16), pmm(16), csx(16), snx(16), ssm(16*nf), asm(16*nf) )
    
    do concurrent ( i2 = 1:16 )
      csx(i2) = roots(i2)
      snx(i2) = sqrt(1-csx(i2)**2)
    end do
    
    do m = 0, this%jmax
      call this%forward_rcb_16_sub( nf, weights(1), sumN(16*nf*m+1), sumS(16*nf*m+1), ssm(1), asm(1) )
      
      !j = m
        mj = m*(this%jmax+1)-(m-2)*(m+1)/2
        
        call this%mmset_16_sub( m, snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_16_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
      
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        
        call this%recursion_16_sub( mj-1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_16_sub( nf, pmj0(1), asm(1), cr(1,mj-1) )
        
        call this%recursion_16_sub( mj, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_16_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
      end do
      
      if ( mod(this%jmax-m,2) /= 0 ) then
        call this%recursion_16_sub( mj+1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_16_sub( nf, pmj0(1), asm(1), cr(1,mj+1) )
      end if
    end do
    
    deallocate( pmj2, pmj1, pmj0, csx, snx, asm, ssm )
    
  end subroutine forward_legesum_16_sub
  
end submodule step16f