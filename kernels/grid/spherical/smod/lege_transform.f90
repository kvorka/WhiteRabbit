submodule (SphericalHarmonics) lege_transform
  implicit none; contains
  
  module pure subroutine lege_transform_sub(this, nforw, nback, cc, cr, grid_sub)
    class(T_lateralGrid),    intent(in)    :: this
    integer,                 intent(in)    :: nforw, nback
    complex(kind=dbl),       intent(in)    :: cc(nback,*)
    complex(kind=dbl),       intent(inout) :: cr(nforw,*)
    integer                                :: i1, i2, i, j, m, mj
    real(kind=dbl),            allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), weight(:), grid(:)
    complex(kind=dbl), pointer             :: pssym(:,:), pasym(:,:), psumN(:,:,:), psumS(:,:,:)
    complex(kind=dbl), target, allocatable :: ssym(:), asym(:), sumN(:), sumS(:)
    
    interface
      module pure subroutine grid_sub(nfour, gxyz)
        integer,                intent(in)    :: nfour
        real(kind=dbl), target, intent(inout) :: gxyz(*)
      end subroutine grid_sub
    end interface
    
    !Allocating needed memory :: step is set to 4
    allocate( pmm(4), pmj(4), pmj1(4), pmj2(4), cosx(4), sinx(4), weight(4), ssym(4*nback), asym(4*nback), &
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
      pssym(1:4,1:nback) => ssym(1:4*nback)
      pasym(1:4,1:nback) => asym(1:4*nback)
      
      psumN(1:nback,1:4,0:this%jmax2) => sumN(1:4*nback*this%jmax3)
      psumS(1:nback,1:4,0:this%jmax2) => sumS(1:4*nback*this%jmax3)
      
      do m = 0, this%jmax2
        call zero_carray_sub( 4*nback, ssym(1) )
        call zero_carray_sub( 4*nback, asym(1) )
        
        if ( m == 0 ) then
          pmm(1:4) = this%cmm(0)
        else
          pmm(1:4) = this%cmm(m) * sinx(1:4) * pmm(1:4)
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:4) = zero
          pmj1(1:4) = zero
          pmj(1:4)  = pmm(1:4)
          
          do concurrent ( i1 = 1:nback, i2 = 1:4 )
            pssym(i2,i1) = pssym(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amj(mj-1) * cosx(i2) * pmj1(i2) - this%bmj(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:4 )
            pasym(i2,i1) = pasym(i2,i1) + cc(i1,mj-1) * pmj(i2)
          end do
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amj(mj) * cosx(i2) * pmj1(i2) - this%bmj(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:4 )
            pssym(i2,i1) = pssym(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amj(mj+1) * cosx(i2) * pmj1(i2) - this%bmj(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:4 )
            pasym(i2,i1) = pasym(i2,i1) + cc(i1,mj+1) * pmj(i2)
          end do
        end if
        
        do concurrent ( i2=1:4, i1=1:nback )
          psumN(i1,i2,m) = pssym(i2,i1) + pasym(i2,i1)
          psumS(i1,i2,m) = pssym(i2,i1) - pasym(i2,i1)
        end do
      end do
      
      do concurrent ( i2=1:4, i1=1:nback )
        psumN(i1,i2,0)%im = zero
        psumS(i1,i2,0)%im = zero
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
      pssym(1:4,1:nforw) => ssym(1:4*nforw)
      pasym(1:4,1:nforw) => asym(1:4*nforw)
      
      psumN(1:nforw,1:4,0:this%jmax2) => sumN(1:4*nforw*this%jmax3)
      psumS(1:nforw,1:4,0:this%jmax2) => sumS(1:4*nforw*this%jmax3)
      
      do m = 0, this%jmax2
        do concurrent ( i2=1:4, i1=1:nforw )
          pssym(i2,i1) = weight(i2) * ( psumN(i1,i2,m) + psumS(i1,i2,m) )
          pasym(i2,i1) = weight(i2) * ( psumN(i1,i2,m) - psumS(i1,i2,m) )
        end do
        
        if ( m == 0 ) then
          do concurrent ( i1=1:nforw, i2=1:4 )
            pssym(i2,i1)%im = zero
            pasym(i2,i1)%im = zero
          end do
          
          pmm(1:4) = this%cmm(0)
        else
          pmm(1:4) = this%cmm(m) * sinx(1:4) * pmm(1:4)
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:4) = zero
          pmj1(1:4) = zero
          pmj(1:4)  = pmm(1:4)
          
          do concurrent ( i1=1:nforw , i2=1:4 )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * pssym(i2,i1)
          end do
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amj(mj-1) * cosx(i2) * pmj1(i2) - this%bmj(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:4 )
            cr(i1,mj-1) = cr(i1,mj-1) + pmj(i2) * pasym(i2,i1)
          end do
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2) = this%amj(mj) * cosx(i2) * pmj1(i2) - this%bmj(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:4 )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * pssym(i2,i1)
          end do
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2) = this%amj(mj+1) * cosx(i2) * pmj1(i2) - this%bmj(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:4 )
            cr(i1,mj+1) = cr(i1,mj+1) + pmj(i2) * pasym(i2,i1)
          end do
        end if
      end do
    end do
    
  end subroutine lege_transform_sub
  
end submodule lege_transform