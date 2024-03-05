submodule(SphericalHarmonics) partial_sums
  implicit none; contains
  
  module pure subroutine partial_backward_2_sub(this, i, n, cosx, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: i, n
    real(kind=dbl),       intent(inout) :: cosx(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(2,*), asym(2,*), sumN(n,2,0:*), sumS(n,2,0:*)
    integer                             :: i1, i2, j, m, mj
    
    cosx(1:2) = this%roots(i:i+1)
    
    call zero_carray_sub( 2*n*this%jmax3, sumN(1,1,0) )
    call zero_carray_sub( 2*n*this%jmax3, sumS(1,1,0) )
    
    do m = 0, this%jmax2
      if ( maxval(abs(this%pmm(i:i+1,m))) < this%tolm ) then
        exit
        
      else
        call zero_carray_sub( 2*n, ssym(1,1) )
        call zero_carray_sub( 2*n, asym(1,1) )
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_backward_set_2_sub( i, m, n, pmj2(1), pmj1(1), pmj(1), cc(1,mj), ssym(1,1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_backward_rec_2_sub( mj-1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj-1), asym(1,1) )
          call this%pmj_backward_rec_2_sub( mj  , n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj  ), ssym(1,1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          call this%pmj_backward_rec_2_sub( mj+1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj+1), asym(1,1) )
        end if
        
        do concurrent ( i2=1:2, i1=1:n )
          sumN(i1,i2,m) = ssym(i2,i1) + asym(i2,i1)
          sumS(i1,i2,m) = ssym(i2,i1) - asym(i2,i1)
        end do
      end if
    end do
    
    call zero_carray_imagpart_sub( 2*n, sumN(1,1,0) )
    call zero_carray_imagpart_sub( 2*n, sumS(1,1,0) )
    
  end subroutine partial_backward_2_sub
  
  module pure subroutine partial_forward_2_sub(this, i, n, weight, cosx, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: i, n
    real(kind=dbl),       intent(in)    :: cosx(*)
    real(kind=dbl),       intent(inout) :: weight(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(2,*), asym(2,*), sumN(n,2,0:*), sumS(n,2,0:*)
    integer                             :: i1, i2, j, m, mj
    
    weight(1:2) = this%fftLege(i:i+1)
    
    do m = 0, this%jmax2
      do concurrent ( i2=1:2, i1=1:n )
        ssym(i2,i1) = weight(i2) * ( sumN(i1,i2,m) + sumS(i1,i2,m) )
        asym(i2,i1) = weight(i2) * ( sumN(i1,i2,m) - sumS(i1,i2,m) )
      end do
      
      if ( m == 0 ) then
        call zero_carray_imagpart_sub( 2*n, ssym(1,1) )
        call zero_carray_imagpart_sub( 2*n, asym(1,1) )
      end if
      
      j = m
        mj = m*this%jmax3-m*(m+1)/2+j+1
        
        call this%pmj_forward_set_2_sub( i, m, n, pmj2(1), pmj1(1), pmj(1), ssym(1,1), cr(1,mj) )
      
      do j = 1, (this%jmax2-m)/2
        mj = mj+2
        
        call this%pmj_forward_rec_2_sub( mj-1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1,1), cr(1,mj-1) )
        call this%pmj_forward_rec_2_sub( mj  , n, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1,1), cr(1,mj  ) )
      end do
      
      if (mod(this%jmax2-m,2) /= 0) then
        call this%pmj_forward_rec_2_sub( mj+1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1,1), cr(1,mj+1) )
      end if
    end do
    
  end subroutine partial_forward_2_sub
  
  module pure subroutine partial_backward_4_sub(this, i, n, cosx, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: i, n
    real(kind=dbl),       intent(inout) :: cosx(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(4,*), asym(4,*), sumN(n,4,0:*), sumS(n,4,0:*)
    integer                             :: i1, i2, j, m, mj
    
    call zero_carray_sub( 4*n*this%jmax3, sumN(1,1,0) )
    call zero_carray_sub( 4*n*this%jmax3, sumS(1,1,0) )
    
    cosx(1:4) = this%roots(i:i+3)
    
    do m = 0, this%jmax2
      if ( maxval(abs(this%pmm(i:i+3,m))) < this%tolm ) then
        exit
        
      else
        call zero_carray_sub( 4*n, ssym(1,1) )
        call zero_carray_sub( 4*n, asym(1,1) )
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_backward_set_4_sub( i, m, n, pmj2(1), pmj1(1), pmj(1), cc(1,mj), ssym(1,1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_backward_rec_4_sub( mj-1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj-1), asym(1,1) )
          call this%pmj_backward_rec_4_sub( mj  , n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj  ), ssym(1,1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          call this%pmj_backward_rec_4_sub( mj+1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj+1), asym(1,1) )
        end if
        
        do concurrent ( i2=1:4, i1=1:n )
          sumN(i1,i2,m) = ssym(i2,i1) + asym(i2,i1)
          sumS(i1,i2,m) = ssym(i2,i1) - asym(i2,i1)
        end do
      end if
    end do
    
    call zero_carray_imagpart_sub( 4*n, sumN(1,1,0) )
    call zero_carray_imagpart_sub( 4*n, sumS(1,1,0) )
    
  end subroutine partial_backward_4_sub
  
  module pure subroutine partial_forward_4_sub(this, i, n, weight, cosx, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: i, n
    real(kind=dbl),       intent(in)    :: cosx(*)
    real(kind=dbl),       intent(inout) :: weight(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(4,*), asym(4,*), sumN(n,4,0:*), sumS(n,4,0:*)
    integer                             :: i1, i2, j, m, mj
    
    weight(1:4) = this%fftLege(i:i+3)
    
    do m = 0, this%jmax2
      do concurrent ( i2=1:4, i1=1:n )
        ssym(i2,i1) = weight(i2) * ( sumN(i1,i2,m) + sumS(i1,i2,m) )
        asym(i2,i1) = weight(i2) * ( sumN(i1,i2,m) - sumS(i1,i2,m) )
      end do
      
      if ( m == 0 ) then
        call zero_carray_imagpart_sub( 4*n, ssym(1,1) )
        call zero_carray_imagpart_sub( 4*n, asym(1,1) )
      end if
      
      j = m
        mj = m*this%jmax3-m*(m+1)/2+j+1
        
        call this%pmj_forward_set_4_sub( i, m, n, pmj2(1), pmj1(1), pmj(1), ssym(1,1), cr(1,mj) )
      
      do j = 1, (this%jmax2-m)/2
        mj = mj+2
        
        call this%pmj_forward_rec_4_sub( mj-1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1,1), cr(1,mj-1) )
        call this%pmj_forward_rec_4_sub( mj  , n, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1,1), cr(1,mj  ) )
      end do
      
      if ( mod(this%jmax2-m,2) /= 0 ) then
        call this%pmj_forward_rec_4_sub( mj+1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1,1), cr(1,mj+1) )
      end if
    end do
    
  end subroutine partial_forward_4_sub
  
  module pure subroutine partial_backward_8_sub(this, i, n, cosx, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: i, n
    real(kind=dbl),       intent(inout) :: cosx(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(8,*), asym(8,*), sumN(n,8,0:*), sumS(n,8,0:*)
    integer                             :: i1, i2, j, m, mj
    
    call zero_carray_sub( 8*n*this%jmax3, sumN(1,1,0) )
    call zero_carray_sub( 8*n*this%jmax3, sumS(1,1,0) )
    
    cosx(1:8) = this%roots(i:i+7)
    
    do m = 0, this%jmax2
      if ( maxval(abs(this%pmm(i:i+7,m))) < this%tolm ) then
        exit
        
      else
        call zero_carray_sub( 8*n, ssym(1,1) )
        call zero_carray_sub( 8*n, asym(1,1) )
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_backward_set_8_sub( i, m, n, pmj2(1), pmj1(1), pmj(1), cc(1,mj), ssym(1,1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_backward_rec_8_sub( mj-1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj-1), asym(1,1) )
          call this%pmj_backward_rec_8_sub( mj  , n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj  ), ssym(1,1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          call this%pmj_backward_rec_8_sub( mj+1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj+1), asym(1,1) )
        end if
        
        do concurrent ( i2=1:8, i1=1:n )
          sumN(i1,i2,m) = ssym(i2,i1) + asym(i2,i1)
          sumS(i1,i2,m) = ssym(i2,i1) - asym(i2,i1)
        end do
      end if
    end do
    
    call zero_carray_imagpart_sub( 8*n, sumN(1,1,0) )
    call zero_carray_imagpart_sub( 8*n, sumS(1,1,0) )
    
  end subroutine partial_backward_8_sub
  
  module pure subroutine partial_forward_8_sub(this, i, n, weight, cosx, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: i, n
    real(kind=dbl),       intent(in)    :: cosx(*)
    real(kind=dbl),       intent(inout) :: weight(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(8,*), asym(8,*), sumN(n,8,0:*), sumS(n,8,0:*)
    integer                             :: i1, i2, j, m, mj
    
    weight(1:8) = this%fftLege(i:i+7)
    
    do m = 0, this%jmax2
      do concurrent ( i2=1:8, i1=1:n )
        ssym(i2,i1) = weight(i2) * ( sumN(i1,i2,m) + sumS(i1,i2,m) )
        asym(i2,i1) = weight(i2) * ( sumN(i1,i2,m) - sumS(i1,i2,m) )
      end do
      
      if ( m == 0 ) then
        call zero_carray_imagpart_sub( 8*n, ssym(1,1) )
        call zero_carray_imagpart_sub( 8*n, asym(1,1) )
      end if
      
      j = m
        mj = m*this%jmax3-m*(m+1)/2+j+1
        
        call this%pmj_forward_set_8_sub( i, m, n, pmj2(1), pmj1(1), pmj(1), ssym(1,1), cr(1,mj) )
      
      do j = 1, (this%jmax2-m)/2
        mj = mj+2
        
        call this%pmj_forward_rec_8_sub( mj-1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1,1), cr(1,mj-1) )
        call this%pmj_forward_rec_8_sub( mj  , n, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1,1), cr(1,mj  ) )
      end do
      
      if ( mod(this%jmax2-m,2) /= 0 ) then
        call this%pmj_forward_rec_8_sub( mj+1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1,1), cr(1,mj+1) )
      end if
    end do
    
  end subroutine partial_forward_8_sub
  
  module pure subroutine partial_backward_16_sub(this, i, n, cosx, pmj2, pmj1, pmj, ssym, asym, cc, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: i, n
    real(kind=dbl),       intent(inout) :: cosx(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(in)    :: cc(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(16,*), asym(16,*), sumN(n,16,0:*), sumS(n,16,0:*)
    integer                             :: i1, i2, j, m, mj
    
    call zero_carray_sub( 16*n*this%jmax3, sumN(1,1,0) )
    call zero_carray_sub( 16*n*this%jmax3, sumS(1,1,0) )
    
    cosx(1:16) = this%roots(i:i+15)
    
    do m = 0, this%jmax2
      if (maxval(abs(this%pmm(i:i+15,m))) < this%tolm) then
        exit
        
      else
        call zero_carray_sub( 16*n, ssym(1,1) )
        call zero_carray_sub( 16*n, asym(1,1) )
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_backward_set_16_sub( i, m, n, pmj2(1), pmj1(1), pmj(1), cc(1,mj), ssym(1,1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_backward_rec_16_sub( mj-1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj-1), asym(1,1) )
          call this%pmj_backward_rec_16_sub( mj  , n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj  ), ssym(1,1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          call this%pmj_backward_rec_16_sub( mj+1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj+1), asym(1,1) )
        end if
        
        do concurrent ( i2=1:16, i1=1:n )
          sumN(i1,i2,m) = ssym(i2,i1) + asym(i2,i1)
          sumS(i1,i2,m) = ssym(i2,i1) - asym(i2,i1)
        end do
      end if
    end do
    
    call zero_carray_imagpart_sub( 16*n, sumN(1,1,0) )
    call zero_carray_imagpart_sub( 16*n, sumS(1,1,0) )
    
  end subroutine partial_backward_16_sub
  
  module pure subroutine partial_forward_16_sub(this, i, n, weight, cosx, pmj2, pmj1, pmj, ssym, asym, cr, sumN, sumS)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: i, n
    real(kind=dbl),       intent(in)    :: cosx(*)
    real(kind=dbl),       intent(inout) :: weight(*), pmj2(*), pmj1(*), pmj(*)
    complex(kind=dbl),    intent(inout) :: cr(n,*)
    complex(kind=dbl),    intent(inout) :: ssym(16,*), asym(16,*), sumN(n,16,0:*), sumS(n,16,0:*)
    integer                             :: i1, i2, j, m, mj
    
    weight(1:16) = this%fftLege(i:i+15)
    
    do m = 0, this%jmax2
      do concurrent ( i2=1:16, i1=1:n )
        ssym(i2,i1) = weight(i2) * ( sumN(i1,i2,m) + sumS(i1,i2,m) )
        asym(i2,i1) = weight(i2) * ( sumN(i1,i2,m) - sumS(i1,i2,m) )
      end do
      
      if ( m == 0 ) then
        call zero_carray_imagpart_sub( 16*n, ssym(1,1) )
        call zero_carray_imagpart_sub( 16*n, asym(1,1) )
      end if
      
      j = m
        mj = m*this%jmax3-m*(m+1)/2+j+1
        
        call this%pmj_forward_set_16_sub( i, m, n, pmj2(1), pmj1(1), pmj(1), ssym(1,1), cr(1,mj) )
      
      do j = 1, (this%jmax2-m)/2
        mj = mj+2
        
        call this%pmj_forward_rec_16_sub( mj-1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1,1), cr(1,mj-1) )
        call this%pmj_forward_rec_16_sub( mj  , n, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1,1), cr(1,mj  ) )
      end do
      
      if ( mod(this%jmax2-m,2) /= 0 ) then
        call this%pmj_forward_rec_16_sub( mj+1, n, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1,1), cr(1,mj+1) )
      end if
    end do
    
  end subroutine partial_forward_16_sub
  
end submodule partial_sums