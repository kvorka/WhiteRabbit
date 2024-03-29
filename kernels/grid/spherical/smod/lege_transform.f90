submodule (SphericalHarmonics) lege_transform
  implicit none; contains
  
  pure subroutine pmj_mmset_sub(m, cmm, snx, pmm, pmj2, pmj1, pmj0)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: cmm, snx(*)
    real(kind=dbl), intent(inout) :: pmm(*)
    real(kind=dbl), intent(out)   :: pmj2(*), pmj1(*), pmj0(*)
    integer                       :: i2
    
    if ( m == 0 ) then
      do concurrent ( i2 = 1:4 )
        pmm(i2) = cmm
      end do
    else
      do concurrent ( i2 = 1:4 )
        pmm(i2) = cmm * snx(i2) * pmm(i2)
      end do
    end if
    
    do concurrent ( i2 = 1:4 )
      pmj2(i2) = zero
      pmj1(i2) = zero
      pmj0(i2) = pmm(i2)
    end do
    
  end subroutine pmj_mmset_sub
  
  pure subroutine pmj_recursion_sub(amj, bmj, csx, pmj2, pmj1, pmj0)
    real(kind=dbl), intent(in)    :: amj, bmj, csx(*)
    real(kind=dbl), intent(inout) :: pmj2(*), pmj1(*), pmj0(*)
    integer                       :: i2
    
    do concurrent ( i2=1:4 )
      pmj2(i2) = pmj1(i2)
      pmj1(i2) = pmj0(i2)
      pmj0(i2) = amj * csx(i2) * pmj1(i2) - bmj * pmj2(i2)
    end do
    
  end subroutine pmj_recursion_sub
  
  pure subroutine pmj_backward_sub(nb, legep, cc, legesum)
    integer,           intent(in)    :: nb
    real(kind=dbl),    intent(in)    :: legep(*)
    complex(kind=dbl), intent(in)    :: cc(*)
    complex(kind=dbl), intent(inout) :: legesum(4,*)
    integer                          :: i1
    complex(kind=dbl)                :: c1
    
    do concurrent ( i1 = 1:nb )
      c1 = cc(i1)
      
      legesum(1,i1) = legesum(1,i1) + c1 * legep(1)
      legesum(2,i1) = legesum(2,i1) + c1 * legep(2)
      legesum(3,i1) = legesum(3,i1) + c1 * legep(3)
      legesum(4,i1) = legesum(4,i1) + c1 * legep(4)
    end do
    
  end subroutine pmj_backward_sub
  
  pure subroutine pmj_forward_sub(nf, legep, legesum, cr)
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: legep(*)
    complex(kind=dbl), intent(in)    :: legesum(4,*)
    complex(kind=dbl), intent(inout) :: cr(*)
    integer                          :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + legep(1) * legesum(1,i1) + legep(2) * legesum(2,i1) + &
                      & legep(3) * legesum(3,i1) + legep(4) * legesum(4,i1)
    end do
    
  end subroutine pmj_forward_sub
  
  module pure subroutine lege_transform_sub(this, nf, nb, cc, cr, grid_sub)
    class(T_lateralGrid),    intent(in)    :: this
    integer,                 intent(in)    :: nf, nb
    complex(kind=dbl),       intent(in)    :: cc(nb,*)
    complex(kind=dbl),       intent(inout) :: cr(nf,*)
    integer                                :: i1, i2, i, j, m, mj
    real(kind=dbl),            allocatable :: pmm(:), pmj0(:), pmj1(:), pmj2(:), csx(:), snx(:), wght(:), grid(:)
    complex(kind=dbl),         allocatable :: ssym(:,:), asym(:,:)
    complex(kind=dbl), pointer             :: psumN(:,:,:), psumS(:,:,:)
    complex(kind=dbl), target, allocatable :: sumN(:), sumS(:)
    
    interface
      module pure subroutine grid_sub(nfour, gxyz)
        integer,                intent(in)    :: nfour
        real(kind=dbl), target, intent(inout) :: gxyz(*)
      end subroutine grid_sub
    end interface
    
    !Allocating memory :: step is set to 4
    allocate( pmm(4), pmj0(4), pmj1(4), pmj2(4), csx(4), snx(4), wght(4), ssym(4,nb), asym(4,nb), &
            & sumN(4*nb*this%jmax3), sumS(4*nb*this%jmax3), grid(4*nb*this%nFourier) )
    
    !Cycle over latitudes :: calculating 4 at once
    do i = 1, this%nLegendre, 4
      do concurrent ( i2 = 1:4 )
        csx(i2)  = this%cosx(i+i2-1)
        snx(i2)  = sqrt(1-csx(i2)**2)
        wght(i2) = this%weight(i+i2-1)
      end do
      
      call zero_carray_sub( 4*nb*this%jmax3, sumN(1) )
      call zero_carray_sub( 4*nb*this%jmax3, sumS(1) )
      
      !Sum of associated Legendre polynomials :: towards grid
      psumN(1:nb,1:4,0:this%jmax2) => sumN(1:4*nb*this%jmax3)
      psumS(1:nb,1:4,0:this%jmax2) => sumS(1:4*nb*this%jmax3)
      
      do m = 0, this%jmax2
        call zero_carray_sub( 4*nb, ssym(1,1) )
        call zero_carray_sub( 4*nb, asym(1,1) )
        
        j = m
          mj = m*this%jmax3-(m-2)*(m+1)/2
          
          call pmj_mmset_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sub( nb, pmj0(1), cc(1,mj), ssym(1,1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call pmj_recursion_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sub( nb, pmj0(1), cc(1,mj-1), asym(1,1) )
          
          call pmj_recursion_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sub( nb, pmj0(1), cc(1,mj), ssym(1,1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          call pmj_recursion_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sub( nb, pmj0(1), cc(1,mj+1), asym(1,1) )
        end if
        
        do concurrent ( i2=1:4, i1=1:nb )
          psumN(i1,i2,m) = ssym(i2,i1) + asym(i2,i1)
          psumS(i1,i2,m) = ssym(i2,i1) - asym(i2,i1)
        end do
      end do
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 4*nb, sumN(1), grid(1) )
      call grid_sub( this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 4*nf, grid(1), sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 4*nb, sumS(1), grid(1) )
      call grid_sub( this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 4*nf, grid(1), sumS(1) )
      
      !Sum of associated Legendre polynomials :: towards space
      psumN(1:nf,1:4,0:this%jmax2) => sumN(1:4*nf*this%jmax3)
      psumS(1:nf,1:4,0:this%jmax2) => sumS(1:4*nf*this%jmax3)
      
      do m = 0, this%jmax2
        do concurrent ( i1=1:nf, i2=1:4 )
          ssym(i2,i1) = wght(i2) * ( psumN(i1,i2,m) + psumS(i1,i2,m) )
          asym(i2,i1) = wght(i2) * ( psumN(i1,i2,m) - psumS(i1,i2,m) )
        end do
        
        j = m
          mj = m*this%jmax3-(m-2)*(m+1)/2
          
          call pmj_mmset_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sub( nf, pmj0(1), ssym(1,1), cr(1,mj) )
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call pmj_recursion_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sub( nf, pmj0(1), asym(1,1), cr(1,mj-1) )
          
          call pmj_recursion_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sub( nf, pmj0(1), ssym(1,1), cr(1,mj) )
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          call pmj_recursion_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sub( nf, pmj0(1), asym(1,1), cr(1,mj+1) )
        end if
      end do
    end do
    
    !Cleaning
    deallocate( pmm, pmj0, pmj1, pmj2, csx, snx, wght, ssym, asym, sumN, sumS, grid )
    
  end subroutine lege_transform_sub
  
end submodule lege_transform