submodule (SphericalHarmonics) lege_transform
  implicit none; contains
  
  pure subroutine pmj_mmset_sub(m, cmm, sinx, pmm, pmj2, pmj1, pmj0)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: cmm, sinx(*)
    real(kind=dbl), intent(inout) :: pmm(*)
    real(kind=dbl), intent(out)   :: pmj2(*), pmj1(*), pmj0(*)
    
    if ( m == 0 ) then
      pmm(1:4) = cmm
    else
      pmm(1:4) = cmm * sinx(1:4) * pmm(1:4)
    end if
    
    pmj2(1:4) = zero
    pmj1(1:4) = zero
    pmj0(1:4) = pmm(1:4)
    
  end subroutine pmj_mmset_sub
  
  pure subroutine pmj_recursion_sub(amj, bmj, cosx, pmj2, pmj1, pmj0)
    real(kind=dbl), intent(in)    :: amj, bmj, cosx(*)
    real(kind=dbl), intent(inout) :: pmj2(*), pmj1(*), pmj0(*)
    integer                       :: i2
    
    do concurrent ( i2=1:4 )
      pmj2(i2) = pmj1(i2)
      pmj1(i2) = pmj0(i2)
      pmj0(i2) = amj * cosx(i2) * pmj1(i2) - bmj * pmj2(i2)
    end do
    
  end subroutine pmj_recursion_sub
  
  pure subroutine pmj_backward_sub(nback, pmj0, cc, legesum)
    integer,           intent(in)    :: nback
    real(kind=dbl),    intent(in)    :: pmj0(*)
    complex(kind=dbl), intent(in)    :: cc(*)
    complex(kind=dbl), intent(inout) :: legesum(4,*)
    integer                          :: i1, i2
    
    do concurrent ( i1 = 1:nback, i2 = 1:4 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * pmj0(i2)
    end do
    
  end subroutine pmj_backward_sub
  
  pure subroutine pmj_forward_sub(nforw, pmj0, legesum, cr)
    integer,           intent(in)    :: nforw
    real(kind=dbl),    intent(in)    :: pmj0(*)
    complex(kind=dbl), intent(in)    :: legesum(4,*)
    complex(kind=dbl), intent(inout) :: cr(*)
    integer                          :: i1
    
    do concurrent ( i1=1:nforw )
      cr(i1) = cr(i1) + sum( pmj0(1:4) * legesum(1:4,i1) )
    end do
    
  end subroutine pmj_forward_sub
  
  module pure subroutine lege_transform_sub(this, nforw, nback, cc, cr, grid_sub)
    class(T_lateralGrid),    intent(in)    :: this
    integer,                 intent(in)    :: nforw, nback
    complex(kind=dbl),       intent(in)    :: cc(nback,*)
    complex(kind=dbl),       intent(inout) :: cr(nforw,*)
    integer                                :: i1, i2, i, j, m, mj
    real(kind=dbl),            allocatable :: pmm(:), pmj0(:), pmj1(:), pmj2(:), cosx(:), sinx(:), weight(:), grid(:)
    complex(kind=dbl),         allocatable :: ssym(:,:), asym(:,:)
    complex(kind=dbl), pointer             :: psumN(:,:,:), psumS(:,:,:)
    complex(kind=dbl), target, allocatable :: sumN(:), sumS(:)
    
    interface
      module pure subroutine grid_sub(nfour, gxyz)
        integer,                intent(in)    :: nfour
        real(kind=dbl), target, intent(inout) :: gxyz(*)
      end subroutine grid_sub
    end interface
    
    !Allocating needed memory :: step is set to 4
    allocate( pmm(4), pmj0(4), pmj1(4), pmj2(4), cosx(4), sinx(4), weight(4), ssym(4,nback), asym(4,nback), &
            & sumN(4*nback*this%jmax3), sumS(4*nback*this%jmax3), grid(4*nback*this%nFourier) )
    
    !Stepping of the algorithm :: 4
    do i = 1, this%nLegendre, 4
      cosx(1:4)   = this%cosx(i:i+3)
      sinx(1:4)   = sqrt(1-cosx(1:4)**2)
      weight(1:4) = this%weight(i:i+3)
      
      call zero_carray_sub( 4*nback*this%jmax3, sumN(1) )
      call zero_carray_sub( 4*nback*this%jmax3, sumS(1) )
      
      !**************************************************************************************************************!
      !The backward (towards grid) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      psumN(1:nback,1:4,0:this%jmax2) => sumN(1:4*nback*this%jmax3)
      psumS(1:nback,1:4,0:this%jmax2) => sumS(1:4*nback*this%jmax3)
      
      do m = 0, this%jmax2
        call zero_carray_sub( 4*nback, ssym(1,1) )
        call zero_carray_sub( 4*nback, asym(1,1) )
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call pmj_mmset_sub( m, this%cmm(m), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sub( nback, pmj0(1), cc(1,mj), ssym(1,1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call pmj_recursion_sub( this%amj(mj-1), this%bmj(mj-1), cosx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sub( nback, pmj0(1), cc(1,mj-1), asym(1,1) )
          
          call pmj_recursion_sub( this%amj(mj), this%bmj(mj), cosx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sub( nback, pmj0(1), cc(1,mj), ssym(1,1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          call pmj_recursion_sub( this%amj(mj+1), this%bmj(mj+1), cosx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sub( nback, pmj0(1), cc(1,mj+1), asym(1,1) )
        end if
        
        do concurrent ( i2=1:4, i1=1:nback )
          psumN(i1,i2,m) = ssym(i2,i1) + asym(i2,i1)
          psumS(i1,i2,m) = ssym(i2,i1) - asym(i2,i1)
        end do
      end do
      
      !**************************************************************************************************************!
      !The backward (towards grid) fft, grid operations and the forward fft (towards space) *************************!
      !**************************************************************************************************************!
      call this%fourtrans%exec_c2r_sub( 4*nback, sumN(1), grid(1) )
      call grid_sub( this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 4*nforw, grid(1), sumN(1) )
      
      call this%fourtrans%exec_c2r_sub( 4*nback, sumS(1), grid(1) )
      call grid_sub( this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 4*nforw, grid(1), sumS(1) )
      
      !**************************************************************************************************************!
      !The forward (towards space) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      psumN(1:nforw,1:4,0:this%jmax2) => sumN(1:4*nforw*this%jmax3)
      psumS(1:nforw,1:4,0:this%jmax2) => sumS(1:4*nforw*this%jmax3)
      
      do m = 0, this%jmax2
        do concurrent ( i2=1:4, i1=1:nforw )
          ssym(i2,i1) = weight(i2) * ( psumN(i1,i2,m) + psumS(i1,i2,m) )
          asym(i2,i1) = weight(i2) * ( psumN(i1,i2,m) - psumS(i1,i2,m) )
        end do
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          call pmj_mmset_sub( m, this%cmm(m), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sub( nforw, pmj0(1), ssym(1,1), cr(1,mj) )
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call pmj_recursion_sub( this%amj(mj-1), this%bmj(mj-1), cosx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sub( nforw, pmj0(1), asym(1,1), cr(1,mj-1) )
          
          call pmj_recursion_sub( this%amj(mj), this%bmj(mj), cosx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sub( nforw, pmj0(1), ssym(1,1), cr(1,mj) )
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          call pmj_recursion_sub( this%amj(mj+1), this%bmj(mj+1), cosx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sub( nforw, pmj0(1), asym(1,1), cr(1,mj+1) )
        end if
      end do
    end do
    
  end subroutine lege_transform_sub
  
end submodule lege_transform