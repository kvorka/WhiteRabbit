submodule (lege_poly) step4f
  implicit none; contains
  
  module pure subroutine forward_sum_4_sub(this, nf, legep, legesum, cr)
    class(T_legep),    intent(in)    :: this
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: legep(4)
    complex(kind=dbl), intent(in)    :: legesum(4,nf)
    complex(kind=dbl), intent(inout) :: cr(nf)
    integer                          :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + sum( legep(1:4) * legesum(1:4,i1) )
    end do
    
  end subroutine forward_sum_4_sub
  
  module pure subroutine forward_rcb_4_sub(this, nf, w, sumN, sumS, sumsym, sumasym)
    class(T_legep),    intent(in)  :: this
    integer,           intent(in)  :: nf
    real(kind=dbl),    intent(in)  :: w(4)
    complex(kind=dbl), intent(in)  :: sumN(nf,4), sumS(nf,4)
    complex(kind=dbl), intent(out) :: sumsym(4,nf), sumasym(4,nf)
    integer                        :: i1, i2
    
    sumsym  = transpose( sumN + sumS )
    sumasym = transpose( sumN - sumS )
    
    do concurrent ( i1 = 1:nf, i2 = 1:4 )
      sumsym(i2,i1)  = w(i2) * sumsym(i2,i1)
      sumasym(i2,i1) = w(i2) * sumasym(i2,i1)
    end do
    
  end subroutine forward_rcb_4_sub
  
  module pure subroutine forward_legesum_4_sub(this, it, nf, sumN, sumS, cr)
    class(T_legep),    intent(in)    :: this
    integer,           intent(in)    :: it, nf
    complex(kind=dbl), intent(in)    :: sumN(*), sumS(*)
    complex(kind=dbl), intent(inout) :: cr(nf,*)
    integer                          :: j, m, mj, i2
    real(kind=dbl),    allocatable   :: pmj2(:), pmj1(:), pmj0(:), pmm(:), csx(:), snx(:), wghts(:)
    complex(kind=dbl), allocatable   :: ssm(:), asm(:)
    
    allocate( pmj2(4), pmj1(4), pmj0(4), pmm(4), csx(4), snx(4), wghts(4), ssm(4*nf), asm(4*nf) )
    
    do concurrent ( i2 = 0:3 )
      csx(i2+1)   = this%roots(it+i2)
      snx(i2+1)   = sqrt(1-this%roots(it+i2)**2)
      wghts(i2+1) = this%weights(it+i2)
    end do
    
    do m = 0, this%jmax
      call this%forward_rcb_4_sub( nf, wghts(1), sumN(4*nf*m+1), sumS(4*nf*m+1), ssm(1), asm(1) )
      
      !j = m
        mj = m*(this%jmax+1)-(m-2)*(m+1)/2
        
        call this%mmset_4_sub( m, snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_4_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
      
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        
        call this%recursion_4_sub( mj-1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_4_sub( nf, pmj0(1), asm(1), cr(1,mj-1) )
        
        call this%recursion_4_sub( mj, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_4_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
      end do
      
      if ( mod(this%jmax-m,2) /= 0 ) then
        call this%recursion_4_sub( mj+1, csx(1), pmj2(1), pmj1(1), pmj0(1) )
        call this%forward_sum_4_sub( nf, pmj0(1), asm(1), cr(1,mj+1) )
      end if
    end do
    
    deallocate( pmj2, pmj1, pmj0, csx, snx, wghts, asm, ssm )
    
  end subroutine forward_legesum_4_sub
  
end submodule step4f