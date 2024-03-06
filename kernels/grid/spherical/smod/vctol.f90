submodule (SphericalHarmonics) vctol
  implicit none ; contains
  
  module subroutine vctol_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer                             :: i, j, m, mj
    real(kind=dbl)                      :: diff, maxdiff, xrand
    real(kind=dbl),       allocatable   :: pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:)
    complex(kind=dbl),    allocatable   :: cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    allocate( cc(this%jms2), cr(this%jms2) )
      
      do m = 0, this%jmax2
        do j = m, this%jmax2
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          if ( m /= 0 ) then
            call random_number( xrand ) ; cc(mj)%re = xrand
            call random_number( xrand ) ; cc(mj)%im = xrand
          else
            call random_number( xrand ) ; cc(mj)%re = xrand
            cc(mj)%im = czero
          end if
        end do
      end do
    
    allocate( pmj(16), pmj1(16), pmj2(16), cosx(16), weight(16), ssym(16), asym(16), &
            & sumN(0:16*this%jmax3-1), sumS(0:16*this%jmax3-1)                       )
      
    do
      cr = czero
      
      !Stepping of the algorithm :: 16
      do i = 1, (this%nLegendre/16)*16, 16
        cosx(1:16) = this%roots(i:i+15)
        
        call zero_carray_sub( 16*this%jmax3, sumN(0) )
        call zero_carray_sub( 16*this%jmax3, sumS(0) )
        
        do m = 0, this%get_maxm_fn(i,16)
          call zero_carray_sub( 16, ssym(1) )
          call zero_carray_sub( 16, asym(1) )
          
          j = m
            mj = m*this%jmax3-m*(m+1)/2+j+1
            
            call this%pmj_backward_set_16_sub( i, m, 1, pmj2(1), pmj1(1), pmj(1), cc(mj), ssym(1) )
            
          do j = 1, (this%jmax2-m)/2
            mj = mj+2
            
            call this%pmj_backward_rec_16_sub( mj-1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj-1), asym(1) )
            call this%pmj_backward_rec_16_sub( mj  , 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj  ), ssym(1) )
          end do
          
          if ( mod((this%jmax2-m),2) /= 0 ) then
            call this%pmj_backward_rec_16_sub( mj+1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj+1), asym(1) )
          end if
          
          call this%pmj_backward_recomb_16_sub( 1, ssym(1), asym(1), sumN(16*m), sumS(16*m) )
        end do
        
        call zero_carray_imagpart_sub( 16, sumN(0) )
        call zero_carray_imagpart_sub( 16, sumS(0) )
        
        weight(1:16) = this%fftLege(i:i+15)
        
        do m = 0, this%jmax2
          call this%pmj_forward_recomb_16_sub( 1, weight(1), sumN(16*m), sumS(16*m), ssym(1), asym(1) )
          
          if ( m == 0 ) then
            call zero_carray_imagpart_sub( 16, ssym(1) )
            call zero_carray_imagpart_sub( 16, asym(1) )
          end if
          
          j = m
            mj = m*this%jmax3-m*(m+1)/2+j+1
            
            call this%pmj_forward_set_16_sub( i, m, 1, pmj2(1), pmj1(1), pmj(1), ssym(1), cr(mj) )
          
          do j = 1, (this%jmax2-m)/2
            mj = mj+2
            
            call this%pmj_forward_rec_16_sub( mj-1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(mj-1) )
            call this%pmj_forward_rec_16_sub( mj  , 1, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1), cr(mj  ) )
          end do
          
          if ( mod(this%jmax2-m,2) /= 0 ) then
            call this%pmj_forward_rec_16_sub( mj+1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(mj+1) )
          end if
        end do
      end do
      
      !Stepping of the algorithm :: 8
      do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
        cosx(1:8) = this%roots(i:i+7)
        
        call zero_carray_sub( 8*this%jmax3, sumN(0) )
        call zero_carray_sub( 8*this%jmax3, sumS(0) )
        
        do m = 0, this%get_maxm_fn(i,8)
          call zero_carray_sub( 8, ssym(1) )
          call zero_carray_sub( 8, asym(1) )
          
          j = m
            mj = m*this%jmax3-m*(m+1)/2+j+1
            
            call this%pmj_backward_set_8_sub( i, m, 1, pmj2(1), pmj1(1), pmj(1), cc(mj), ssym(1) )
            
          do j = 1, (this%jmax2-m)/2
            mj = mj+2
            
            call this%pmj_backward_rec_8_sub( mj-1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj-1), asym(1) )
            call this%pmj_backward_rec_8_sub( mj  , 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj  ), ssym(1) )
          end do
          
          if ( mod((this%jmax2-m),2) /= 0 ) then
            call this%pmj_backward_rec_8_sub( mj+1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj+1), asym(1) )
          end if
          
          call this%pmj_backward_recomb_8_sub( 1, ssym(1), asym(1), sumN(8*m), sumS(8*m) )
        end do
        
        call zero_carray_imagpart_sub( 8, sumN(0) )
        call zero_carray_imagpart_sub( 8, sumS(0) )
        
        weight(1:8) = this%fftLege(i:i+7)
        
        do m = 0, this%jmax2
          call this%pmj_forward_recomb_8_sub( 1, weight(1), sumN(8*m), sumS(8*m), ssym(1), asym(1) )
          
          if ( m == 0 ) then
            call zero_carray_imagpart_sub( 8, ssym(1) )
            call zero_carray_imagpart_sub( 8, asym(1) )
          end if
          
          j = m
            mj = m*this%jmax3-m*(m+1)/2+j+1
            
            call this%pmj_forward_set_8_sub( i, m, 1, pmj2(1), pmj1(1), pmj(1), ssym(1), cr(mj) )
          
          do j = 1, (this%jmax2-m)/2
            mj = mj+2
            
            call this%pmj_forward_rec_8_sub( mj-1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(mj-1) )
            call this%pmj_forward_rec_8_sub( mj  , 1, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1), cr(mj  ) )
          end do
          
          if ( mod(this%jmax2-m,2) /= 0 ) then
            call this%pmj_forward_rec_8_sub( mj+1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(mj+1) )
          end if
        end do
      end do
      
      !Stepping of the algorithm :: 4
      do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
        cosx(1:4) = this%roots(i:i+3)
        
        call zero_carray_sub( 4*this%jmax3, sumN(0) )
        call zero_carray_sub( 4*this%jmax3, sumS(0) )
        
        do m = 0, this%get_maxm_fn(i,4)
          call zero_carray_sub( 4, ssym(1) )
          call zero_carray_sub( 4, asym(1) )
          
          j = m
            mj = m*this%jmax3-m*(m+1)/2+j+1
            
            call this%pmj_backward_set_4_sub( i, m, 1, pmj2(1), pmj1(1), pmj(1), cc(mj), ssym(1) )
            
          do j = 1, (this%jmax2-m)/2
            mj = mj+2
            
            call this%pmj_backward_rec_4_sub( mj-1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj-1), asym(1) )
            call this%pmj_backward_rec_4_sub( mj  , 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj  ), ssym(1) )
          end do
          
          if ( mod((this%jmax2-m),2) /= 0 ) then
            call this%pmj_backward_rec_4_sub( mj+1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj+1), asym(1) )
          end if
          
          call this%pmj_backward_recomb_4_sub( 1, ssym(1), asym(1), sumN(4*m), sumS(4*m) )
        end do
        
        call zero_carray_imagpart_sub( 4, sumN(0) )
        call zero_carray_imagpart_sub( 4, sumS(0) )
        
        weight(1:4) = this%fftLege(i:i+3)
        
        do m = 0, this%jmax2
          call this%pmj_forward_recomb_4_sub( 1, weight(1), sumN(4*m), sumS(4*m), ssym(1), asym(1) )
          
          if ( m == 0 ) then
            call zero_carray_imagpart_sub( 4, ssym(1) )
            call zero_carray_imagpart_sub( 4, asym(1) )
          end if
          
          j = m
            mj = m*this%jmax3-m*(m+1)/2+j+1
            
            call this%pmj_forward_set_4_sub( i, m, 1, pmj2(1), pmj1(1), pmj(1), ssym(1), cr(mj) )
          
          do j = 1, (this%jmax2-m)/2
            mj = mj+2
            
            call this%pmj_forward_rec_4_sub( mj-1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(mj-1) )
            call this%pmj_forward_rec_4_sub( mj  , 1, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1), cr(mj  ) )
          end do
          
          if ( mod(this%jmax2-m,2) /= 0 ) then
            call this%pmj_forward_rec_4_sub( mj+1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(mj+1) )
          end if
        end do
      end do
      
      !Stepping of the algorithm :: 2
      do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
        cosx(1:2) = this%roots(i:i+1)
        
        call zero_carray_sub( 2*this%jmax3, sumN(0) )
        call zero_carray_sub( 2*this%jmax3, sumS(0) )
        
        do m = 0, this%get_maxm_fn(i,2)
          call zero_carray_sub( 2, ssym(1) )
          call zero_carray_sub( 2, asym(1) )
          
          j = m
            mj = m*this%jmax3-m*(m+1)/2+j+1
            
            call this%pmj_backward_set_2_sub( i, m, 1, pmj2(1), pmj1(1), pmj(1), cc(mj), ssym(1) )
            
          do j = 1, (this%jmax2-m)/2
            mj = mj+2
            
            call this%pmj_backward_rec_2_sub( mj-1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj-1), asym(1) )
            call this%pmj_backward_rec_2_sub( mj  , 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj  ), ssym(1) )
          end do
          
          if ( mod((this%jmax2-m),2) /= 0 ) then
            call this%pmj_backward_rec_2_sub( mj+1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(mj+1), asym(1) )
          end if
          
          call this%pmj_backward_recomb_2_sub( 1, ssym(1), asym(1), sumN(2*m), sumS(2*m) )
        end do
        
        call zero_carray_imagpart_sub( 2, sumN(0) )
        call zero_carray_imagpart_sub( 2, sumS(0) )
        
        weight(1:2) = this%fftLege(i:i+1)
        
        do m = 0, this%jmax2
          call this%pmj_forward_recomb_2_sub( 1, weight(1), sumN(2*m), sumS(2*m), ssym(1), asym(1) )
          
          if ( m == 0 ) then
            call zero_carray_imagpart_sub( 2, ssym(1) )
            call zero_carray_imagpart_sub( 2, asym(1) )
          end if
          
          j = m
            mj = m*this%jmax3-m*(m+1)/2+j+1
            
            call this%pmj_forward_set_2_sub( i, m, 1, pmj2(1), pmj1(1), pmj(1), ssym(1), cr(mj) )
          
          do j = 1, (this%jmax2-m)/2
            mj = mj+2
            
            call this%pmj_forward_rec_2_sub( mj-1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(mj-1) )
            call this%pmj_forward_rec_2_sub( mj  , 1, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1), cr(mj  ) )
          end do
          
          if ( mod(this%jmax2-m,2) /= 0 ) then
            call this%pmj_forward_rec_2_sub( mj+1, 1, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(mj+1) )
          end if
        end do
      end do
        
      maxdiff = 0._dbl
      
      do mj = 1, this%jms2
        diff = abs( abs( cc(mj) / cr(mj) ) / this%scale / sqrt(4*pi) - 1 )
          
        if ( diff > maxdiff ) maxdiff = diff
      end do
      
      if ( this%tolm < 1.0d-100) then
        this%tolm = 1.0d-90
        exit
      else if ( maxdiff <= 1.0d-4 ) then
        exit
      else
        this%tolm = this%tolm / 10
      end if
    end do
      
    deallocate( cr, cc, sumN, sumS, pmj, pmj1, pmj2, cosx, weight, ssym, asym )
    
  end subroutine vctol_sub
  
end submodule vctol