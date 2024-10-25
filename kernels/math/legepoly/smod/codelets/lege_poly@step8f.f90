submodule (lege_poly) step8f
  implicit none; contains
  
  module pure subroutine forward_sum_8_sub(this, nf, legep, legesum, cr)
    class(T_legep),    intent(in)    :: this
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: legep(8)
    complex(kind=dbl), intent(in)    :: legesum(8,nf)
    complex(kind=dbl), intent(inout) :: cr(nf)
    integer                          :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + sum( legep(1:8) * legesum(1:8,i1) )
    end do
    
  end subroutine forward_sum_8_sub
  
  module pure subroutine forward_rcb_8_sub(this, nf, w, sumN, sumS, sumsym, sumasym)
    class(T_legep),    intent(in)  :: this
    integer,           intent(in)  :: nf
    real(kind=dbl),    intent(in)  :: w(8)
    real(kind=dbl),    intent(in)  :: sumN(8,nf,2), sumS(8,nf,2)
    complex(kind=dbl), intent(out) :: sumsym(8,nf), sumasym(8,nf)
    integer                        :: i1, i2
    
    do concurrent ( i2 = 1:nf, i1 = 1:8 )
      sumsym(i1,i2)  = cmplx( sumN(i1,i2,1) + sumS(i1,i2,1), sumN(i1,i2,2) + sumS(i1,i2,2), kind=dbl )
      sumasym(i1,i2) = cmplx( sumN(i1,i2,1) - sumS(i1,i2,1), sumN(i1,i2,2) - sumS(i1,i2,2), kind=dbl )
    end do
    
    do concurrent ( i2 = 1:nf, i1 = 1:8 )
      sumsym(i1,i2)  = w(i1) * sumsym(i1,i2)
      sumasym(i1,i2) = w(i1) * sumasym(i1,i2)
    end do
    
  end subroutine forward_rcb_8_sub
  
  module pure subroutine forward_legesum_8_sub(this, it, nf, sumN, sumS, cr)
    class(T_legep),    intent(in)    :: this
    integer,           intent(in)    :: it, nf
    real(kind=dbl),    intent(in)    :: sumN(8,nf,2,0:this%jmax), sumS(8,nf,2,0:this%jmax)
    complex(kind=dbl), intent(inout) :: cr(nf,*)
    integer                          :: j, m, mj, i2
    real(kind=dbl),    allocatable   :: pmj2(:), pmj1(:), pmj0(:), pmm(:), csx(:), snx(:), wgx(:)
    complex(kind=dbl), allocatable   :: ssm(:), asm(:)
    
    allocate( pmj2(8), pmj1(8), pmj0(8), pmm(8), csx(8), snx(8), wgx(8), ssm(8*nf), asm(8*nf) )
    
    do concurrent ( i2 = 0:7 )
      csx(i2+1) = this%rootsweights(1,it+i2)
      snx(i2+1) = this%rootsweights(2,it+i2)
      wgx(i2+1) = this%rootsweights(3,it+i2)
    end do
    
    do m = 0, this%jmax
      call this%forward_rcb_8_sub( nf, wgx(1), sumN(1,1,1,m), sumS(1,1,1,m), ssm(1), asm(1) )
      
      !j = m
        mj = m*(this%jmax+1)-(m-2)*(m+1)/2
        
        call this%mmset_8_sub( m, snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_8_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
      
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        
        call this%recursion_8_sub( mj-1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_8_sub( nf, pmj0(1), asm(1), cr(1,mj-1) )
        
        call this%recursion_8_sub( mj, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_8_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
      end do
      
      if ( mod(this%jmax-m,2) /= 0 ) then
        call this%recursion_8_sub( mj+1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_8_sub( nf, pmj0(1), asm(1), cr(1,mj+1) )
      end if
    end do
    
    deallocate( pmj2, pmj1, pmj0, pmm, csx, snx, wgx, asm, ssm )
    
  end subroutine forward_legesum_8_sub
  
end submodule step8f