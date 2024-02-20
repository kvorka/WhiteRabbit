submodule(SphericalHarmonics) partial_sums
  implicit none; contains
  
  module pure subroutine partial_backward_2_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(n,2,*), sumS(n,2,*)
    integer                             :: j, m, mj
        
    do m = 0, this%maxj
      call pmm_recursion_2_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm(1:2))) < this%tolm) exit
      
      call zero_carray_sub( 2*n, ssym(1) )
      call zero_carray_sub( 2*n, asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j+1
        
        call pmj_setup_2_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_2_sub( n, pmj(1), cc(1,mj), ssym(1) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_2_sub( this%amjrr(mj-1), this%bmjrr(mj-1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_2_sub( n, pmj(1), cc(1,mj-1), asym(1) )
        
        call pmj_recursion_2_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_2_sub( n, pmj(1), cc(1,mj), ssym(1) )
        
        if ( maxval(abs(pmj(1:2))) < this%tolm ) exit
      end do
      
      if ( (maxval(abs(pmj(1:2))) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
        mj = mj+1
        
        call pmj_recursion_2_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_2_sub( n, pmj(1), cc(1,mj), asym(1) )
      end if
      
      call pmj_backward_recomb_2_sub( m, n, ssym(1), asym(1), sumN(1,1,m+1), sumS(1,1,m+1) )
    end do
    
  end subroutine partial_backward_2_sub
  
  module pure subroutine partial_forward_2_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(n,2,*), sumS(n,2,*)
    integer                             :: j, m, mj
    
    do m = 0, this%maxj
      call pmm_recursion_2_sub( m, sinx(1), pmm(1) )
      call pmj_forward_recomb_2_sub( m, n, weight(1), sumN(1,1,m+1), sumS(1,1,m+1), ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j+1
        
        call pmj_setup_2_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_2_sub( n, pmj(1), ssym(1), cr(1,mj) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_2_sub( this%amjrr(mj-1), this%bmjrr(mj-1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_2_sub( n, pmj(1), asym(1), cr(1,mj-1) )
        
        call pmj_recursion_2_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_2_sub( n, pmj(1), ssym(1), cr(1,mj) )
      end do
      
      if (mod(this%maxj-m,2) /= 0) then
        mj = mj+1
        
        call pmj_recursion_2_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_2_sub( n, pmj(1), asym(1), cr(1,mj) )
      end if
    end do
    
  end subroutine partial_forward_2_sub
  
  module pure subroutine partial_backward_4_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(n,4,*), sumS(n,4,*)
    integer                             :: j, m, mj
        
    do m = 0, this%maxj
      call pmm_recursion_4_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm(1:4))) < this%tolm) exit
      
      call zero_carray_sub( 4*n, ssym(1) )
      call zero_carray_sub( 4*n, asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j+1
        
        call pmj_setup_4_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_4_sub( n, pmj(1), cc(1,mj), ssym(1) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_4_sub( this%amjrr(mj-1), this%bmjrr(mj-1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_4_sub( n, pmj(1), cc(1,mj-1), asym(1) )
        
        call pmj_recursion_4_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_4_sub( n, pmj(1), cc(1,mj), ssym(1) )
        
        if ( maxval(abs(pmj(1:4))) < this%tolm ) exit
      end do
      
      if ( (maxval(abs(pmj(1:4))) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
        mj = mj+1
        
        call pmj_recursion_4_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_4_sub( n, pmj(1), cc(1,mj), asym(1) )
      end if
      
      call pmj_backward_recomb_4_sub( m, n, ssym(1), asym(1), sumN(1,1,m+1), sumS(1,1,m+1) )
    end do
    
  end subroutine partial_backward_4_sub
  
  module pure subroutine partial_forward_4_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(n,4,*), sumS(n,4,*)
    integer                             :: j, m, mj
    
    do m = 0, this%maxj
      call pmm_recursion_4_sub( m, sinx(1), pmm(1) )
      call pmj_forward_recomb_4_sub( m, n, weight(1), sumN(1,1,m+1), sumS(1,1,m+1), ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j+1
        
        call pmj_setup_4_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_4_sub( n, pmj(1), ssym(1), cr(1,mj) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_4_sub( this%amjrr(mj-1), this%bmjrr(mj-1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_4_sub( n, pmj(1), asym(1), cr(1,mj-1) )
        
        call pmj_recursion_4_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_4_sub( n, pmj(1), ssym(1), cr(1,mj) )
      end do
      
      if (mod(this%maxj-m,2) /= 0) then
        mj = mj+1
        
        call pmj_recursion_4_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_4_sub( n, pmj(1), asym(1), cr(1,mj) )
      end if
    end do
    
  end subroutine partial_forward_4_sub
  
  module pure subroutine partial_backward_8_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(n,8,*), sumS(n,8,*)
    integer                             :: j, m, mj
    
    do m = 0, this%maxj
      call pmm_recursion_8_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm(1:8))) < this%tolm) exit
      
      call zero_carray_sub( 8*n, ssym(1) )
      call zero_carray_sub( 8*n, asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j+1
        
        call pmj_setup_8_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_8_sub( n, pmj(1), cc(1,mj), ssym(1) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_8_sub( this%amjrr(mj-1), this%bmjrr(mj-1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_8_sub( n, pmj(1), cc(1,mj-1), asym(1) )
        
        call pmj_recursion_8_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_8_sub( n, pmj(1), cc(1,mj), ssym(1) )
        
        if ( maxval(abs(pmj(1:8))) < this%tolm ) exit
      end do
      
      if ( (maxval(abs(pmj(1:8))) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
        mj = mj+1
        
        call pmj_recursion_8_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_8_sub( n, pmj(1), cc(1,mj), asym(1) )
      end if
      
      call pmj_backward_recomb_8_sub( m, n, ssym(1), asym(1), sumN(1,1,m+1), sumS(1,1,m+1) )
    end do
    
  end subroutine partial_backward_8_sub
  
  module pure subroutine partial_forward_8_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(n,8,*), sumS(n,8,*)
    integer                             :: j, m, mj
    
    do m = 0, this%maxj
      call pmm_recursion_8_sub( m, sinx(1), pmm(1) )
      call pmj_forward_recomb_8_sub( m, n, weight(1), sumN(1,1,m+1), sumS(1,1,m+1), ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j+1
        
        call pmj_setup_8_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_8_sub( n, pmj(1), ssym(1), cr(1,mj) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_8_sub( this%amjrr(mj-1), this%bmjrr(mj-1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_8_sub( n, pmj(1), asym(1), cr(1,mj-1) )
        
        call pmj_recursion_8_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_8_sub( n, pmj(1), ssym(1), cr(1,mj) )
      end do
      
      if (mod(this%maxj-m,2) /= 0) then
        mj = mj+1
        
        call pmj_recursion_8_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_8_sub( n, pmj(1), asym(1), cr(1,mj) )
      end if
    end do
    
  end subroutine partial_forward_8_sub
  
  module pure subroutine partial_backward_16_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(n,16,*), sumS(n,16,*)
    integer                             :: j, m, mj
        
    do m = 0, this%maxj
      call pmm_recursion_16_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm(1:16))) < this%tolm) exit
      
      call zero_carray_sub( 16*n, ssym(1) )
      call zero_carray_sub( 16*n, asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j+1
        
        call pmj_setup_16_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_16_sub( n, pmj(1), cc(1,mj), ssym(1) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_16_sub( this%amjrr(mj-1), this%bmjrr(mj-1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_16_sub( n, pmj(1), cc(1,mj-1), asym(1) )
        
        call pmj_recursion_16_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_16_sub( n, pmj(1), cc(1,mj), ssym(1) )
        
        if ( maxval(abs(pmj(1:16))) < this%tolm ) exit
      end do
      
      if ( (maxval(abs(pmj(1:16))) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
        mj = mj+1
        
        call pmj_recursion_16_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_16_sub( n, pmj(1), cc(1,mj), asym(1) )
      end if
      
      call pmj_backward_recomb_16_sub( m, n, ssym(1), asym(1), sumN(1,1,m+1), sumS(1,1,m+1) )
    end do
    
  end subroutine partial_backward_16_sub
  
  module pure subroutine partial_forward_16_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(n,16,*), sumS(n,16,*)
    integer                             :: j, m, mj
    
    do m = 0, this%maxj
      call pmm_recursion_16_sub( m, sinx(1), pmm(1) )
      call pmj_forward_recomb_16_sub( m, n, weight(1), sumN(1,1,m+1), sumS(1,1,m+1), ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j+1
        
        call pmj_setup_16_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_16_sub( n, pmj(1), ssym(1), cr(1,mj) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_16_sub( this%amjrr(mj-1), this%bmjrr(mj-1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_16_sub( n, pmj(1), asym(1), cr(1,mj-1) )
        
        call pmj_recursion_16_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_16_sub( n, pmj(1), ssym(1), cr(1,mj) )
      end do
      
      if (mod(this%maxj-m,2) /= 0) then
        mj = mj+1
        
        call pmj_recursion_16_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_16_sub( n, pmj(1), asym(1), cr(1,mj) )
      end if
    end do
    
  end subroutine partial_forward_16_sub
  
end submodule partial_sums