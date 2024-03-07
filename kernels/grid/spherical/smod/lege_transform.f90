submodule (SphericalHarmonics) lege_transform
  implicit none; contains
  
  module pure subroutine lege_transform_sub(this, nforw, nback, cc, cr, grid_2_sub, grid_4_sub, grid_8_sub, grid_16_sub)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: nforw, nback
    complex(kind=dbl),    intent(in)    :: cc(nback,*)
    complex(kind=dbl),    intent(inout) :: cr(nforw,*)
    integer                             :: i, j, m, mj
    real(kind=dbl),       allocatable   :: pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:), grid(:)
    complex(kind=dbl),    allocatable   :: ssym(:), asym(:), sumN(:), sumS(:)
    
    interface
      pure subroutine grid_2_sub(sph, gxyz, sumNS)
        import dbl, T_lateralGrid
        class(T_lateralGrid),   intent(in)    :: sph
        real(kind=dbl), target, intent(out)   :: gxyz(*)
        complex(kind=dbl),      intent(inout) :: sumNS(*)
      end subroutine grid_2_sub
      
      pure subroutine grid_4_sub(sph, gxyz, sumNS)
        import dbl, T_lateralGrid
        class(T_lateralGrid),   intent(in)    :: sph
        real(kind=dbl), target, intent(out)   :: gxyz(*)
        complex(kind=dbl),      intent(inout) :: sumNS(*)
      end subroutine grid_4_sub
      
      pure subroutine grid_8_sub(sph, gxyz, sumNS)
        import dbl, T_lateralGrid
        class(T_lateralGrid),   intent(in)    :: sph
        real(kind=dbl), target, intent(out)   :: gxyz(*)
        complex(kind=dbl),      intent(inout) :: sumNS(*)
      end subroutine grid_8_sub
      
      pure subroutine grid_16_sub(sph, gxyz, sumNS)
        import dbl, T_lateralGrid
        class(T_lateralGrid),   intent(in)    :: sph
        real(kind=dbl), target, intent(out)   :: gxyz(*)
        complex(kind=dbl),      intent(inout) :: sumNS(*)
      end subroutine grid_16_sub
    end interface
    
    !Allocating needed memory :: no reallocate for lower stepping
    allocate( pmj(16), pmj1(16), pmj2(16), cosx(16), weight(16), ssym(16*nback), asym(16*nback),         &
            & sumN(0:16*nback*this%jmax3-1), sumS(0:16*nback*this%jmax3-1), grid(16*nback*this%nFourier) )
    
    !Stepping of the algorithm :: 16
    do i = 1, (this%nLegendre/16)*16, 16
      cosx(1:16) = this%roots(i:i+15)
      
      call zero_carray_sub( 16*nback*this%jmax3, sumN(0) )
      call zero_carray_sub( 16*nback*this%jmax3, sumS(0) )
      
      do m = 0, this%get_maxm_fn(i,16)
        call zero_carray_sub( 16*nback, ssym(1) )
        call zero_carray_sub( 16*nback, asym(1) )
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_backward_set_16_sub( i, m, nback, pmj2(1), pmj1(1), pmj(1), cc(1,mj), ssym(1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_backward_rec_16_sub( mj-1, nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj-1), asym(1) )
          call this%pmj_backward_rec_16_sub( mj  , nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj  ), ssym(1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          call this%pmj_backward_rec_16_sub( mj+1, nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj+1), asym(1) )
        end if
        
        call this%pmj_backward_recomb_16_sub( nback, ssym(1), asym(1), sumN(16*nback*m), sumS(16*nback*m) )
      end do
      
      call zero_carray_imagpart_sub( 16*nback, sumN(0) )
      call zero_carray_imagpart_sub( 16*nback, sumS(0) )
      
      call grid_16_sub( this, grid(1), sumN(0) )
      call grid_16_sub( this, grid(1), sumS(0) )
      
      weight(1:16) = this%fftLege(i:i+15)
      
      do m = 0, this%jmax2
        call this%pmj_forward_recomb_16_sub( nforw, weight(1), sumN(16*nforw*m), sumS(16*nforw*m), ssym(1), asym(1) )
        
        if ( m == 0 ) then
          call zero_carray_imagpart_sub( 16*nforw, ssym(1) )
          call zero_carray_imagpart_sub( 16*nforw, asym(1) )
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_forward_set_16_sub( i, m, nforw, pmj2(1), pmj1(1), pmj(1), ssym(1), cr(1,mj) )
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_forward_rec_16_sub( mj-1, nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(1,mj-1) )
          call this%pmj_forward_rec_16_sub( mj  , nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1), cr(1,mj  ) )
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          call this%pmj_forward_rec_16_sub( mj+1, nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(1,mj+1) )
        end if
      end do
    end do
    
    !Stepping of the algorithm :: 8
    do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
      cosx(1:8) = this%roots(i:i+7)
      
      call zero_carray_sub( 8*nback*this%jmax3, sumN(0) )
      call zero_carray_sub( 8*nback*this%jmax3, sumS(0) )
      
      do m = 0, this%get_maxm_fn(i,8)
        call zero_carray_sub( 8*nback, ssym(1) )
        call zero_carray_sub( 8*nback, asym(1) )
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_backward_set_8_sub( i, m, nback, pmj2(1), pmj1(1), pmj(1), cc(1,mj), ssym(1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_backward_rec_8_sub( mj-1, nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj-1), asym(1) )
          call this%pmj_backward_rec_8_sub( mj  , nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj  ), ssym(1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          call this%pmj_backward_rec_8_sub( mj+1, nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj+1), asym(1) )
        end if
        
        call this%pmj_backward_recomb_8_sub( nback, ssym(1), asym(1), sumN(8*nback*m), sumS(8*nback*m) )
      end do
      
      call zero_carray_imagpart_sub( 8*nback, sumN(0) )
      call zero_carray_imagpart_sub( 8*nback, sumS(0) )
      
      call grid_8_sub( this, grid(1), sumN(0) )
      call grid_8_sub( this, grid(1), sumS(0) )
      
      weight(1:8) = this%fftLege(i:i+7)
      
      do m = 0, this%jmax2
        call this%pmj_forward_recomb_8_sub( nforw, weight(1), sumN(8*nforw*m), sumS(8*nforw*m), ssym(1), asym(1) )
        
        if ( m == 0 ) then
          call zero_carray_imagpart_sub( 8*nforw, ssym(1) )
          call zero_carray_imagpart_sub( 8*nforw, asym(1) )
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_forward_set_8_sub( i, m, nforw, pmj2(1), pmj1(1), pmj(1), ssym(1), cr(1,mj) )
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_forward_rec_8_sub( mj-1, nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(1,mj-1) )
          call this%pmj_forward_rec_8_sub( mj  , nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1), cr(1,mj  ) )
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          call this%pmj_forward_rec_8_sub( mj+1, nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(1,mj+1) )
        end if
      end do
    end do
    
    !Stepping of the algorithm :: 4
    do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
      cosx(1:4) = this%roots(i:i+3)
      
      call zero_carray_sub( 4*nback*this%jmax3, sumN(0) )
      call zero_carray_sub( 4*nback*this%jmax3, sumS(0) )
      
      do m = 0, this%get_maxm_fn(i,4)
        call zero_carray_sub( 4*nback, ssym(1) )
        call zero_carray_sub( 4*nback, asym(1) )
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_backward_set_4_sub( i, m, nback, pmj2(1), pmj1(1), pmj(1), cc(1,mj), ssym(1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_backward_rec_4_sub( mj-1, nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj-1), asym(1) )
          call this%pmj_backward_rec_4_sub( mj  , nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj  ), ssym(1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          call this%pmj_backward_rec_4_sub( mj+1, nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj+1), asym(1) )
        end if
        
        call this%pmj_backward_recomb_4_sub( nback, ssym(1), asym(1), sumN(4*nback*m), sumS(4*nback*m) )
      end do
      
      call zero_carray_imagpart_sub( 4*nback, sumN(0) )
      call zero_carray_imagpart_sub( 4*nback, sumS(0) )
      
      call grid_4_sub( this, grid(1), sumN(0) )
      call grid_4_sub( this, grid(1), sumS(0) )
      
      weight(1:4) = this%fftLege(i:i+3)
      
      do m = 0, this%jmax2
        call this%pmj_forward_recomb_4_sub( nforw, weight(1), sumN(4*nforw*m), sumS(4*nforw*m), ssym(1), asym(1) )
        
        if ( m == 0 ) then
          call zero_carray_imagpart_sub( 4*nforw, ssym(1) )
          call zero_carray_imagpart_sub( 4*nforw, asym(1) )
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_forward_set_4_sub( i, m, nforw, pmj2(1), pmj1(1), pmj(1), ssym(1), cr(1,mj) )
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_forward_rec_4_sub( mj-1, nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(1,mj-1) )
          call this%pmj_forward_rec_4_sub( mj  , nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1), cr(1,mj  ) )
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          call this%pmj_forward_rec_4_sub( mj+1, nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(1,mj+1) )
        end if
      end do
    end do
    
    !Stepping of the algorithm :: 2
    do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
      cosx(1:2) = this%roots(i:i+1)
      
      call zero_carray_sub( 2*nback*this%jmax3, sumN(0) )
      call zero_carray_sub( 2*nback*this%jmax3, sumS(0) )
      
      do m = 0, this%get_maxm_fn(i,2)
        call zero_carray_sub( 2*nback, ssym(1) )
        call zero_carray_sub( 2*nback, asym(1) )
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_backward_set_2_sub( i, m, nback, pmj2(1), pmj1(1), pmj(1), cc(1,mj), ssym(1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_backward_rec_2_sub( mj-1, nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj-1), asym(1) )
          call this%pmj_backward_rec_2_sub( mj  , nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj  ), ssym(1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          call this%pmj_backward_rec_2_sub( mj+1, nback, cosx(1), pmj2(1), pmj1(1), pmj(1), cc(1,mj+1), asym(1) )
        end if
        
        call this%pmj_backward_recomb_2_sub( nback, ssym(1), asym(1), sumN(2*nback*m), sumS(2*nback*m) )
      end do
      
      call zero_carray_imagpart_sub( 2*nback, sumN(0) )
      call zero_carray_imagpart_sub( 2*nback, sumS(0) )
      
      call grid_2_sub( this, grid(1), sumN(0) )
      call grid_2_sub( this, grid(1), sumS(0) )
      
      weight(1:2) = this%fftLege(i:i+1)
      
      do m = 0, this%jmax2
        call this%pmj_forward_recomb_2_sub( nforw, weight(1), sumN(2*nforw*m), sumS(2*nforw*m), ssym(1), asym(1) )
        
        if ( m == 0 ) then
          call zero_carray_imagpart_sub( 2*nforw, ssym(1) )
          call zero_carray_imagpart_sub( 2*nforw, asym(1) )
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call this%pmj_forward_set_2_sub( i, m, nforw, pmj2(1), pmj1(1), pmj(1), ssym(1), cr(1,mj) )
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call this%pmj_forward_rec_2_sub( mj-1, nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(1,mj-1) )
          call this%pmj_forward_rec_2_sub( mj  , nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), ssym(1), cr(1,mj  ) )
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          call this%pmj_forward_rec_2_sub( mj+1, nforw, cosx(1), pmj2(1), pmj1(1), pmj(1), asym(1), cr(1,mj+1) )
        end if
      end do
    end do
    
    deallocate( sumN, sumS, grid, pmj, pmj1, pmj2, cosx, weight, ssym, asym )
    
  end subroutine lege_transform_sub
  
end submodule lege_transform