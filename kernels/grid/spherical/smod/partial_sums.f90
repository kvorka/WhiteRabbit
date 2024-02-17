submodule(SphericalHarmonics) partial_sums
  implicit none; contains
  
  module pure subroutine partial_backward_2_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    integer                             :: j, m, mj
        
    do m = 0, this%maxj
      call pmm_recursion_2_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm(1:2))) < this%tolm) exit
      
      call zero_poly_sub( 2*n, ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j
        
        call pmj_setup_2_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_2_sub( n, pmj(1), cc(1+n*mj), ssym(1) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_2_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_2_sub( n, pmj(1), cc(1+n*(mj-1)), asym(1) )
        
        call pmj_recursion_2_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_2_sub( n, pmj(1), cc(1+n*mj), ssym(1) )
        
        if ( maxval(abs(pmj(1:2))) < this%tolm ) exit
      end do
      
      if ( (maxval(abs(pmj(1:2))) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
        mj = mj+1
        
        call pmj_recursion_2_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_2_sub( n, pmj(1), cc(1+n*mj), asym(1) )
      end if
      
      call pmj_backward_recomb_2_sub( m, n, ssym(1), asym(1), sumN(1+2*n*m), sumS(1+2*n*m) )
    end do
    
  end subroutine partial_backward_2_sub
  
  module pure subroutine partial_forward_2_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    integer                             :: j, m, mj
    
    do m = 0, this%maxj
      call pmm_recursion_2_sub( m, sinx(1), pmm(1) )
      call pmj_forward_recomb_2_sub( m, n, weight(1), sumN(1+2*n*m), sumS(1+2*n*m), ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j
        
        call pmj_setup_2_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_2_sub( n, pmj(1), ssym(1), cr(1+n*mj) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_2_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_2_sub( n, pmj(1), asym(1), cr(1+n*(mj-1)) )
        
        call pmj_recursion_2_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_2_sub( n, pmj(1), ssym(1), cr(1+n*mj) )
      end do
      
      if (mod(this%maxj-m,2) /= 0) then
        mj = mj+1
        
        call pmj_recursion_2_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_2_sub( n, pmj(1), asym(1), cr(1+n*mj) )
      end if
    end do
    
  end subroutine partial_forward_2_sub
  
  module pure subroutine partial_backward_4_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    integer                             :: j, m, mj
        
    do m = 0, this%maxj
      call pmm_recursion_4_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm(1:4))) < this%tolm) exit
      
      call zero_poly_sub( 4*n, ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j
        
        call pmj_setup_4_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_4_sub( n, pmj(1), cc(1+n*mj), ssym(1) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_4_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_4_sub( n, pmj(1), cc(1+n*(mj-1)), asym(1) )
        
        call pmj_recursion_4_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_4_sub( n, pmj(1), cc(1+n*mj), ssym(1) )
        
        if ( maxval(abs(pmj(1:4))) < this%tolm ) exit
      end do
      
      if ( (maxval(abs(pmj(1:4))) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
        mj = mj+1
        
        call pmj_recursion_4_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_4_sub( n, pmj(1), cc(1+n*mj), asym(1) )
      end if
      
      call pmj_backward_recomb_4_sub( m, n, ssym(1), asym(1), sumN(1+4*n*m), sumS(1+4*n*m) )
    end do
    
  end subroutine partial_backward_4_sub
  
  module pure subroutine partial_forward_4_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    integer                             :: j, m, mj
    
    do m = 0, this%maxj
      call pmm_recursion_4_sub( m, sinx(1), pmm(1) )
      call pmj_forward_recomb_4_sub( m, n, weight(1), sumN(1+4*n*m), sumS(1+4*n*m), ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j
        
        call pmj_setup_4_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_4_sub( n, pmj(1), ssym(1), cr(1+n*mj) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_4_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_4_sub( n, pmj(1), asym(1), cr(1+n*(mj-1)) )
        
        call pmj_recursion_4_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_4_sub( n, pmj(1), ssym(1), cr(1+n*mj) )
      end do
      
      if (mod(this%maxj-m,2) /= 0) then
        mj = mj+1
        
        call pmj_recursion_4_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_4_sub( n, pmj(1), asym(1), cr(1+n*mj) )
      end if
    end do
    
  end subroutine partial_forward_4_sub
  
  module pure subroutine partial_backward_8_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    integer                             :: j, m, mj
        
    do m = 0, this%maxj
      call pmm_recursion_8_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm(1:8))) < this%tolm) exit
      
      call zero_poly_sub( 8*n, ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j
        
        call pmj_setup_8_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_8_sub( n, pmj(1), cc(1+n*mj), ssym(1) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_8_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_8_sub( n, pmj(1), cc(1+n*(mj-1)), asym(1) )
        
        call pmj_recursion_8_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_8_sub( n, pmj(1), cc(1+n*mj), ssym(1) )
        
        if ( maxval(abs(pmj(1:8))) < this%tolm ) exit
      end do
      
      if ( (maxval(abs(pmj(1:8))) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
        mj = mj+1
        
        call pmj_recursion_8_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_8_sub( n, pmj(1), cc(1+n*mj), asym(1) )
      end if
      
      call pmj_backward_recomb_8_sub( m, n, ssym(1), asym(1), sumN(1+8*n*m), sumS(1+8*n*m) )
    end do
    
  end subroutine partial_backward_8_sub
  
  module pure subroutine partial_forward_8_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    integer                             :: j, m, mj
    
    do m = 0, this%maxj
      call pmm_recursion_8_sub( m, sinx(1), pmm(1) )
      call pmj_forward_recomb_8_sub( m, n, weight(1), sumN(1+8*n*m), sumS(1+8*n*m), ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j
        
        call pmj_setup_8_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_8_sub( n, pmj(1), ssym(1), cr(1+n*mj) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_8_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_8_sub( n, pmj(1), asym(1), cr(1+n*(mj-1)) )
        
        call pmj_recursion_8_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_8_sub( n, pmj(1), ssym(1), cr(1+n*mj) )
      end do
      
      if (mod(this%maxj-m,2) /= 0) then
        mj = mj+1
        
        call pmj_recursion_8_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_8_sub( n, pmj(1), asym(1), cr(1+n*mj) )
      end if
    end do
    
  end subroutine partial_forward_8_sub
  
  module pure subroutine partial_backward_16_sub(this, n, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    integer                             :: j, m, mj
        
    do m = 0, this%maxj
      call pmm_recursion_16_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm(1:16))) < this%tolm) exit
      
      call zero_poly_sub( 16*n, ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j
        
        call pmj_setup_16_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_16_sub( n, pmj(1), cc(1+n*mj), ssym(1) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_16_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_16_sub( n, pmj(1), cc(1+n*(mj-1)), asym(1) )
        
        call pmj_recursion_16_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_16_sub( n, pmj(1), cc(1+n*mj), ssym(1) )
        
        if ( maxval(abs(pmj(1:16))) < this%tolm ) exit
      end do
      
      if ( (maxval(abs(pmj(1:16))) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
        mj = mj+1
        
        call pmj_recursion_16_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_backward_16_sub( n, pmj(1), cc(1+n*mj), asym(1) )
      end if
      
      call pmj_backward_recomb_16_sub( m, n, ssym(1), asym(1), sumN(1+16*n*m), sumS(1+16*n*m) )
    end do
    
  end subroutine partial_backward_16_sub
  
  module pure subroutine partial_forward_16_sub(this, n, weight, cosx, sinx, pmm, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: n
    real(kind=dbl),       intent(in)    :: cosx(*), sinx(*), weight(*)
    real(kind=dbl),       intent(inout) :: pmm(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(*)
    complex(kind=dbl),    intent(inout) :: ssym(*), asym(*), sumN(*), sumS(*)
    integer                             :: j, m, mj
    
    do m = 0, this%maxj
      call pmm_recursion_16_sub( m, sinx(1), pmm(1) )
      call pmj_forward_recomb_16_sub( m, n, weight(1), sumN(1+16*n*m), sumS(1+16*n*m), ssym(1), asym(1) )
      
      j = m
        mj = m*(this%maxj+1)-m*(m+1)/2+j
        
        call pmj_setup_16_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_16_sub( n, pmj(1), ssym(1), cr(1+n*mj) )
      
      do j = 1, (this%maxj-m)/2
        mj = mj+2
        
        call pmj_recursion_16_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_16_sub( n, pmj(1), asym(1), cr(1+n*(mj-1)) )
        
        call pmj_recursion_16_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_16_sub( n, pmj(1), ssym(1), cr(1+n*mj) )
      end do
      
      if (mod(this%maxj-m,2) /= 0) then
        mj = mj+1
        
        call pmj_recursion_16_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
        call pmj_forward_16_sub( n, pmj(1), asym(1), cr(1+n*mj) )
      end if
    end do
    
  end subroutine partial_forward_16_sub
  
end submodule partial_sums