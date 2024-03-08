submodule (SphericalHarmonics) lege_transform
  implicit none; contains
  
  module pure subroutine rescale_sub( this, arr, length )
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: length
    complex(kind=dbl),    intent(inout) :: arr(*)
    integer                             :: ijm
    
    do concurrent ( ijm = 1:length )
      arr(ijm) = this%scale * arr(ijm)
    end do
    
  end subroutine rescale_sub
  
  module pure integer function get_maxm_fn(this, i, i2)
    class(T_lateralGrid), intent(in) :: this
    integer,              intent(in) :: i, i2
    integer                          :: m, maxm
    
    if ( allocated(this%maxm) ) then
      maxm = this%maxm(i)
      
    else
      maxm = this%jmax2
      
      do m = 0, this%jmax2
        if ( maxval(abs(this%pmm(i:i+i2-1,m))) < this%tolm ) then
          maxm = m-1
          exit
        end if
      end do
    end if
    
    get_maxm_fn = maxm
    
  end function get_maxm_fn
  
  module pure subroutine lege_transform_sub(this, nforw, nback, cc, cr, grid_sub)
    class(T_lateralGrid),    intent(in)    :: this
    integer,                 intent(in)    :: nforw, nback
    complex(kind=dbl),       intent(in)    :: cc(nback,*)
    complex(kind=dbl),       intent(inout) :: cr(nforw,*)
    integer                                :: i1, i2, i, j, m, mj
    real(kind=dbl),            allocatable :: pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:), grid(:)
    complex(kind=dbl), pointer             :: pssym(:,:), pasym(:,:), psumN(:,:,:), psumS(:,:,:)
    complex(kind=dbl), target, allocatable :: ssym(:), asym(:), sumN(:), sumS(:)
    
    interface
      pure subroutine grid_sub(nfour, nstep, gxyz)
        import dbl, T_lateralGrid
        integer,                intent(in)    :: nfour, nstep
        real(kind=dbl), target, intent(inout) :: gxyz(*)
      end subroutine grid_sub
    end interface
    
    !Allocating needed memory :: no reallocate for lower stepping
    allocate( pmj(16), pmj1(16), pmj2(16), cosx(16), weight(16), ssym(16*nback), asym(16*nback), &
            & sumN(16*nback*this%jmax3), sumS(16*nback*this%jmax3), grid(16*nback*this%nFourier) )
    
    !Stepping of the algorithm :: 16
    do i = 1, (this%nLegendre/16)*16, 16
      !**************************************************************************************************************!
      !Array preparations *******************************************************************************************!
      !**************************************************************************************************************!
      cosx(1:16)   = this%roots(i:i+15)
      weight(1:16) = this%fftLege(i:i+15)
      
      !**************************************************************************************************************!
      !The backward (towards grid) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      pssym(1:16,1:nback) => ssym(1:16*nback)
      pasym(1:16,1:nback) => asym(1:16*nback)
      
      psumN(1:nback,1:16,0:this%jmax2) => sumN(1:16*nback*this%jmax3)
      psumS(1:nback,1:16,0:this%jmax2) => sumS(1:16*nback*this%jmax3)
      
      do concurrent ( m=0:this%jmax2, i2=1:16, i1=1:nback )
        psumN(i1,i2,m) = czero
        psumS(i1,i2,m) = czero
      end do
      
      do m = 0, this%get_maxm_fn(i,16)
        do concurrent ( i1=1:nback, i2=1:16 )
          pssym(i2,i1) = czero
          pasym(i2,i1) = czero
        end do
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:16) = zero
          pmj1(1:16) = zero
          pmj(1:16)  = this%pmm(i:i+15,m)
          
          do concurrent ( i1 = 1:nback, i2 = 1:16 )
            pssym(i2,i1) = pssym(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:16) = pmj1(1:16)
          pmj1(1:16) = pmj(1:16)
          
          do concurrent ( i2=1:16 )
            pmj(i2)  = this%amjrr(mj-1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:16 )
            pasym(i2,i1) = pasym(i2,i1) + cc(i1,mj-1) * pmj(i2)
          end do
          
          pmj2(1:16) = pmj1(1:16)
          pmj1(1:16) = pmj(1:16)
          
          do concurrent ( i2=1:16 )
            pmj(i2)  = this%amjrr(mj) * cosx(i2) * pmj1(i2) - this%bmjrr(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:16 )
            pssym(i2,i1) = pssym(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          pmj2(1:16) = pmj1(1:16)
          pmj1(1:16) = pmj(1:16)
          
          do concurrent ( i2=1:16 )
            pmj(i2)  = this%amjrr(mj+1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:16 )
            pasym(i2,i1) = pasym(i2,i1) + cc(i1,mj+1) * pmj(i2)
          end do
        end if
        
        do concurrent ( i2=1:16, i1=1:nback )
          psumN(i1,i2,m) = pssym(i2,i1) + pasym(i2,i1)
          psumS(i1,i2,m) = pssym(i2,i1) - pasym(i2,i1)
        end do
      end do
      
      do concurrent ( i2=1:16, i1=1:nback )
        psumN(i1,i2,0)%im = zero
        psumS(i1,i2,0)%im = zero
      end do
      
      !**************************************************************************************************************!
      !The backward (towards grid) fft, grid operations and the forward fft (towards space) *************************!
      !**************************************************************************************************************!
      call this%fourtrans%exec_c2r_sub( 16*nback, sumN(1), grid(1) )
      call grid_sub( this%nFourier, 16, grid(1) )
      call this%fourtrans%exec_r2c_sub( 16*nforw, grid(1), sumN(1) )
      
      call this%fourtrans%exec_c2r_sub( 16*nback, sumS(1), grid(1) )
      call grid_sub( this%nFourier, 16, grid(1) )
      call this%fourtrans%exec_r2c_sub( 16*nforw, grid(1), sumS(1) )
      
      !**************************************************************************************************************!
      !The forward (towards space) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      pssym(1:16,1:nforw) => ssym(1:16*nforw)
      pasym(1:16,1:nforw) => asym(1:16*nforw)
      
      psumN(1:nforw,1:16,0:this%jmax2) => sumN(1:16*nforw*this%jmax3)
      psumS(1:nforw,1:16,0:this%jmax2) => sumS(1:16*nforw*this%jmax3)
      
      do m = 0, this%jmax2
        do concurrent ( i2=1:16, i1=1:nforw )
          pssym(i2,i1) = weight(i2) * ( psumN(i1,i2,m) + psumS(i1,i2,m) )
          pasym(i2,i1) = weight(i2) * ( psumN(i1,i2,m) - psumS(i1,i2,m) )
        end do
        
        if ( m == 0 ) then
          do concurrent ( i1=1:nforw, i2=1:16 )
            pssym(i2,i1)%im = zero
            pasym(i2,i1)%im = zero
          end do
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:16) = zero
          pmj1(1:16) = zero
          pmj(1:16)  = this%pmm(i:i+15,m)
          
          do concurrent ( i1=1:nforw , i2=1:16 )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * pssym(i2,i1)
          end do
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:16) = pmj1(1:16)
          pmj1(1:16) = pmj(1:16)
          
          do concurrent ( i2=1:16 )
            pmj(i2)  = this%amjrr(mj-1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:16 )
            cr(i1,mj-1) = cr(i1,mj-1) + pmj(i2) * pasym(i2,i1)
          end do

          pmj2(1:16) = pmj1(1:16)
          pmj1(1:16) = pmj(1:16)
          
          do concurrent ( i2=1:16 )
            pmj(i2)  = this%amjrr(mj) * cosx(i2) * pmj1(i2) - this%bmjrr(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:16 )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * pssym(i2,i1)
          end do
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          pmj2(1:16) = pmj1(1:16)
          pmj1(1:16) = pmj(1:16)
          
          do concurrent ( i2=1:16 )
            pmj(i2)  = this%amjrr(mj+1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:16 )
            cr(i1,mj+1) = cr(i1,mj+1) + pmj(i2) * pasym(i2,i1)
          end do
        end if
      end do
    end do
    
    !Stepping of the algorithm :: 8
    do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
      !**************************************************************************************************************!
      !Array preparations *******************************************************************************************!
      !**************************************************************************************************************!
      cosx(1:8)   = this%roots(i:i+7)
      weight(1:8) = this%fftLege(i:i+7)
      
      !**************************************************************************************************************!
      !The backward (towards grid) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      pssym(1:8,1:nback) => ssym(1:8*nback)
      pasym(1:8,1:nback) => asym(1:8*nback)
      
      psumN(1:nback,1:8,0:this%jmax2) => sumN(1:8*nback*this%jmax3)
      psumS(1:nback,1:8,0:this%jmax2) => sumS(1:8*nback*this%jmax3)
      
      do concurrent ( m=0:this%jmax2, i2=1:8, i1=1:nback )
        psumN(i1,i2,m) = czero
        psumS(i1,i2,m) = czero
      end do
      
      do m = 0, this%get_maxm_fn(i,8)
        do concurrent ( i1=1:nback, i2=1:8 )
          pssym(i2,i1) = czero
          pasym(i2,i1) = czero
        end do
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:8) = zero
          pmj1(1:8) = zero
          pmj(1:8)  = this%pmm(i:i+7,m)
          
          do concurrent ( i1 = 1:nback, i2 = 1:8 )
            pssym(i2,i1) = pssym(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amjrr(mj-1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:8 )
            pasym(i2,i1) = pasym(i2,i1) + cc(i1,mj-1) * pmj(i2)
          end do
          
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amjrr(mj) * cosx(i2) * pmj1(i2) - this%bmjrr(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:8 )
            pssym(i2,i1) = pssym(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amjrr(mj+1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:8 )
            pasym(i2,i1) = pasym(i2,i1) + cc(i1,mj+1) * pmj(i2)
          end do
        end if
        
        do concurrent ( i2=1:8, i1=1:nback )
          psumN(i1,i2,m) = pssym(i2,i1) + pasym(i2,i1)
          psumS(i1,i2,m) = pssym(i2,i1) - pasym(i2,i1)
        end do
      end do
      
      do concurrent ( i2=1:8, i1=1:nback )
        psumN(i1,i2,0)%im = zero
        psumS(i1,i2,0)%im = zero
      end do
      
      !**************************************************************************************************************!
      !The backward (towards grid) fft, grid operations and the forward fft (towards space) *************************!
      !**************************************************************************************************************!
      call this%fourtrans%exec_c2r_sub( 8*nback, sumN(1), grid(1) )
      call grid_sub( this%nFourier, 8, grid(1) )
      call this%fourtrans%exec_r2c_sub( 8*nforw, grid(1), sumN(1) )
      
      call this%fourtrans%exec_c2r_sub( 8*nback, sumS(1), grid(1) )
      call grid_sub( this%nFourier, 8, grid(1) )
      call this%fourtrans%exec_r2c_sub( 8*nforw, grid(1), sumS(1) )
      
      !**************************************************************************************************************!
      !The forward (towards space) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      pssym(1:8,1:nforw) => ssym(1:8*nforw)
      pasym(1:8,1:nforw) => asym(1:8*nforw)
      
      psumN(1:nforw,1:8,0:this%jmax2) => sumN(1:8*nforw*this%jmax3)
      psumS(1:nforw,1:8,0:this%jmax2) => sumS(1:8*nforw*this%jmax3)
      
      do m = 0, this%jmax2
        do concurrent ( i2=1:8, i1=1:nforw )
          pssym(i2,i1) = weight(i2) * ( psumN(i1,i2,m) + psumS(i1,i2,m) )
          pasym(i2,i1) = weight(i2) * ( psumN(i1,i2,m) - psumS(i1,i2,m) )
        end do
        
        if ( m == 0 ) then
          do concurrent ( i1=1:nforw, i2=1:8 )
            pssym(i2,i1)%im = zero
            pasym(i2,i1)%im = zero
          end do
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:8) = zero
          pmj1(1:8) = zero
          pmj(1:8)  = this%pmm(i:i+7,m)
          
          do concurrent ( i1=1:nforw , i2=1:8 )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * pssym(i2,i1)
          end do
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amjrr(mj-1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:8 )
            cr(i1,mj-1) = cr(i1,mj-1) + pmj(i2) * pasym(i2,i1)
          end do
          
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amjrr(mj) * cosx(i2) * pmj1(i2) - this%bmjrr(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:8 )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * pssym(i2,i1)
          end do
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2) = this%amjrr(mj+1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:8 )
            cr(i1,mj+1) = cr(i1,mj+1) + pmj(i2) * pasym(i2,i1)
          end do
        end if
      end do
    end do
    
    !Stepping of the algorithm :: 4
    do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
      !**************************************************************************************************************!
      !Array preparations *******************************************************************************************!
      !**************************************************************************************************************!
      cosx(1:4)   = this%roots(i:i+3)
      weight(1:4) = this%fftLege(i:i+3)
      
      !**************************************************************************************************************!
      !The backward (towards grid) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      pssym(1:4,1:nback) => ssym(1:4*nback)
      pasym(1:4,1:nback) => asym(1:4*nback)
      
      psumN(1:nback,1:4,0:this%jmax2) => sumN(1:4*nback*this%jmax3)
      psumS(1:nback,1:4,0:this%jmax2) => sumS(1:4*nback*this%jmax3)
      
      do concurrent ( m=0:this%jmax2, i2=1:4, i1=1:nback )
        psumN(i1,i2,m) = czero
        psumS(i1,i2,m) = czero
      end do
      
      do m = 0, this%get_maxm_fn(i,4)
        do concurrent ( i1=1:nback, i2=1:4 )
          pssym(i2,i1) = czero
          pasym(i2,i1) = czero
        end do
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:4) = zero
          pmj1(1:4) = zero
          pmj(1:4)  = this%pmm(i:i+3,m)
          
          do concurrent ( i1 = 1:nback, i2 = 1:4 )
            pssym(i2,i1) = pssym(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amjrr(mj-1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:4 )
            pasym(i2,i1) = pasym(i2,i1) + cc(i1,mj-1) * pmj(i2)
          end do
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amjrr(mj) * cosx(i2) * pmj1(i2) - this%bmjrr(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:4 )
            pssym(i2,i1) = pssym(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amjrr(mj+1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj+1) * pmj2(i2)
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
      call grid_sub( this%nFourier, 4, grid(1) )
      call this%fourtrans%exec_r2c_sub( 4*nforw, grid(1), sumN(1) )
      
      call this%fourtrans%exec_c2r_sub( 4*nback, sumS(1), grid(1) )
      call grid_sub( this%nFourier, 4, grid(1) )
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
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:4) = zero
          pmj1(1:4) = zero
          pmj(1:4)  = this%pmm(i:i+3,m)
          
          do concurrent ( i1=1:nforw , i2=1:4 )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * pssym(i2,i1)
          end do
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amjrr(mj-1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:4 )
            cr(i1,mj-1) = cr(i1,mj-1) + pmj(i2) * pasym(i2,i1)
          end do
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2) = this%amjrr(mj) * cosx(i2) * pmj1(i2) - this%bmjrr(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:4 )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * pssym(i2,i1)
          end do
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2) = this%amjrr(mj+1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:4 )
            cr(i1,mj+1) = cr(i1,mj+1) + pmj(i2) * pasym(i2,i1)
          end do
        end if
      end do
    end do
    
    !Stepping of the algorithm :: 2
    do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
      !**************************************************************************************************************!
      !Array preparations *******************************************************************************************!
      !**************************************************************************************************************!
      cosx(1:2)   = this%roots(i:i+1)
      weight(1:2) = this%fftLege(i:i+1)
      
      !**************************************************************************************************************!
      !The backward (towards grid) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      pssym(1:2,1:nback) => ssym(1:2*nback)
      pasym(1:2,1:nback) => asym(1:2*nback)
      
      psumN(1:nback,1:2,0:this%jmax2) => sumN(1:2*nback*this%jmax3)
      psumS(1:nback,1:2,0:this%jmax2) => sumS(1:2*nback*this%jmax3)
      
      do concurrent ( m=0:this%jmax2, i2=1:2, i1=1:nback )
        psumN(i1,i2,m) = czero
        psumS(i1,i2,m) = czero
      end do
      
      do m = 0, this%get_maxm_fn(i,2)
        do concurrent ( i1=1:nback, i2=1:2 )
          pssym(i2,i1) = czero
          pasym(i2,i1) = czero
        end do
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:2) = zero
          pmj1(1:2) = zero
          pmj(1:2)  = this%pmm(i:i+1,m)
          
          do concurrent ( i1 = 1:nback, i2 = 1:2 )
            pssym(i2,i1) = pssym(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:2) = pmj1(1:2)
          pmj1(1:2) = pmj(1:2)
          
          do concurrent ( i2=1:2 )
            pmj(i2)  = this%amjrr(mj-1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:2 )
            pasym(i2,i1) = pasym(i2,i1) + cc(i1,mj-1) * pmj(i2)
          end do
          
          pmj2(1:2) = pmj1(1:2)
          pmj1(1:2) = pmj(1:2)
          
          do concurrent ( i2=1:2 )
            pmj(i2)  = this%amjrr(mj) * cosx(i2) * pmj1(i2) - this%bmjrr(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:2 )
            pssym(i2,i1) = pssym(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          pmj2(1:2) = pmj1(1:2)
          pmj1(1:2) = pmj(1:2)
          
          do concurrent ( i2=1:2 )
            pmj(i2)  = this%amjrr(mj+1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:2 )
            pasym(i2,i1) = pasym(i2,i1) + cc(i1,mj+1) * pmj(i2)
          end do
        end if
        
        do concurrent ( i2=1:2, i1=1:nback )
          psumN(i1,i2,m) = pssym(i2,i1) + pasym(i2,i1)
          psumS(i1,i2,m) = pssym(i2,i1) - pasym(i2,i1)
        end do
      end do
      
      do concurrent ( i2=1:2, i1=1:nback )
        psumN(i1,i2,0)%im = zero
        psumS(i1,i2,0)%im = zero
      end do
      
      !**************************************************************************************************************!
      !The backward (towards grid) fft, grid operations and the forward fft (towards space) *************************!
      !**************************************************************************************************************!
      call this%fourtrans%exec_c2r_sub( 2*nback, sumN(1), grid(1) )
      call grid_sub( this%nFourier, 2, grid(1) )
      call this%fourtrans%exec_r2c_sub( 2*nforw, grid(1), sumN(1) )
      
      call this%fourtrans%exec_c2r_sub( 2*nback, sumS(1), grid(1) )
      call grid_sub( this%nFourier, 2, grid(1) )
      call this%fourtrans%exec_r2c_sub( 2*nforw, grid(1), sumS(1) )
      
      !**************************************************************************************************************!
      !The forward (towards space) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      pssym(1:2,1:nforw) => ssym(1:2*nforw)
      pasym(1:2,1:nforw) => asym(1:2*nforw)
      
      psumN(1:nforw,1:2,0:this%jmax2) => sumN(1:2*nforw*this%jmax3)
      psumS(1:nforw,1:2,0:this%jmax2) => sumS(1:2*nforw*this%jmax3)
      
      do m = 0, this%jmax2
        do concurrent ( i2=1:2, i1=1:nforw )
          pssym(i2,i1) = weight(i2) * ( psumN(i1,i2,m) + psumS(i1,i2,m) )
          pasym(i2,i1) = weight(i2) * ( psumN(i1,i2,m) - psumS(i1,i2,m) )
        end do
        
        if ( m == 0 ) then
          do concurrent ( i1=1:nforw, i2=1:2 )
            pssym(i2,i1)%im = zero
            pasym(i2,i1)%im = zero
          end do
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:2) = zero
          pmj1(1:2) = zero
          pmj(1:2)  = this%pmm(i:i+1,m)
          
          do concurrent ( i1=1:nforw , i2=1:2 )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * pssym(i2,i1)
          end do
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:2) = pmj1(1:2)
          pmj1(1:2) = pmj(1:2)
          
          do concurrent ( i2=1:2 )
            pmj(i2) = this%amjrr(mj-1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:2 )
            cr(i1,mj-1) = cr(i1,mj-1) + pmj(i2) * pasym(i2,i1)
          end do
          
          pmj2(1:2) = pmj1(1:2)
          pmj1(1:2) = pmj(1:2)
          
          do concurrent ( i2=1:2 )
            pmj(i2)  = this%amjrr(mj) * cosx(i2) * pmj1(i2) - this%bmjrr(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:2 )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * pssym(i2,i1)
          end do
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          pmj2(1:2) = pmj1(1:2)
          pmj1(1:2) = pmj(1:2)
          
          do concurrent ( i2=1:2 )
            pmj(i2)  = this%amjrr(mj+1) * cosx(i2) * pmj1(i2) - this%bmjrr(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:2 )
            cr(i1,mj+1) = cr(i1,mj+1) + pmj(i2) * pasym(i2,i1)
          end do
        end if
      end do
    end do
    
    deallocate( sumN, sumS, grid, pmj, pmj1, pmj2, cosx, weight, ssym, asym )
    
  end subroutine lege_transform_sub
  
end submodule lege_transform