submodule (SphericalHarmonics) lege_transform
  implicit none; contains
  
  module subroutine lege_transform_sub(this, nforw, nback, cc, cr, grid_sub)
    class(T_lateralGrid),      intent(in)    :: this
    integer,                   intent(in)    :: nforw, nback
    complex(kind=dbl), target, intent(in)    :: cc(nback,*)
    complex(kind=dbl), target, intent(inout) :: cr(nforw,*)
    integer                                  :: i1, i2, i3, i, j, m, mj
    real(kind=dbl),            allocatable   :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), weight(:)
    real(kind=dbl), pointer                  :: pcc(:,:,:), pcr(:,:,:), pssym(:,:,:), pasym(:,:,:), psumN(:,:,:,:), psumS(:,:,:,:)
    real(kind=dbl), target,    allocatable   :: ssym(:), asym(:), sumN(:), sumS(:)
    type(c_ptr)                              :: cptr
    
    interface
      module pure subroutine grid_sub(nfour, gxyz)
        integer,                intent(in)    :: nfour
        real(kind=dbl), target, intent(inout) :: gxyz(*)
      end subroutine grid_sub
    end interface
    
    !Preparing pointers to input and output
    cptr = c_loc(cc); call c_f_pointer(cptr, pcc, shape=[2,nback,this%jms2])
    cptr = c_loc(cr); call c_f_pointer(cptr, pcr, shape=[2,nforw,this%jms2])
    
    !Allocating needed memory :: step is set to 8
    allocate( pmm(8), pmj(8), pmj1(8), pmj2(8), cosx(8), sinx(8), weight(8), ssym(16*nback), &
            & asym(16*nback), sumN(8*nback*this%nFourier), sumS(8*nback*this%nFourier)       )
    
    !Stepping of the algorithm :: 8
    do i = 1, this%nLegendre, 8
      cosx(1:8)   = this%cosx(i:i+7)
      sinx(1:8)   = sqrt(1-cosx(1:8)**2)
      weight(1:8) = this%weight(i:i+7)
      
      sumN = zero
      sumS = zero
      
      !**************************************************************************************************************!
      !The backward (towards grid) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      pssym(1:8,1:2,1:nback) => ssym(1:16*nback)
      pasym(1:8,1:2,1:nback) => asym(1:16*nback)
      
      psumN(1:nback,1:8,1:2,0:this%jmax2) => sumN(1:16*nback*this%jmax3)
      psumS(1:nback,1:8,1:2,0:this%jmax2) => sumS(1:16*nback*this%jmax3)
      
      do m = 0, this%jmax2
        ssym = zero
        asym = zero
        
        if ( m == 0 ) then
          pmm(1:8) = this%cmm(0)
        else
          pmm(1:8) = this%cmm(m) * sinx(1:8) * pmm(1:8)
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:8) = zero
          pmj1(1:8) = zero
          pmj(1:8)  = pmm(1:8)
          
          do concurrent ( i1 = 1:nback, i2 = 1:2, i3 = 1:8 )
            pssym(i3,i2,i1) = pssym(i3,i2,i1) + pcc(i2,i1,mj) * pmj(i3)
          end do
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amj(mj-1) * cosx(i2) * pmj1(i2) - this%bmj(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1 = 1:nback, i2 = 1:2, i3=1:8 )
            pasym(i3,i2,i1) = pasym(i3,i2,i1) + pcc(i2,i1,mj-1) * pmj(i3)
          end do
          
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amj(mj) * cosx(i2) * pmj1(i2) - this%bmj(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1 = 1:nback, i2 = 1:2, i3=1:8 )
            pssym(i3,i2,i1) = pssym(i3,i2,i1) + pcc(i2,i1,mj) * pmj(i3)
          end do
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amj(mj+1) * cosx(i2) * pmj1(i2) - this%bmj(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1 = 1:nback, i2 = 1:2, i3=1:8 )
            pasym(i3,i2,i1) = pasym(i3,i2,i1) + pcc(i2,i1,mj+1) * pmj(i3)
          end do
        end if
        
        do concurrent ( i2=1:8, i1=1:nback, i3=1:2 )
          psumN(i1,i2,i3,m) = pssym(i2,i3,i1) + pasym(i2,i3,i1)
          psumS(i1,i2,i3,m) = pssym(i2,i3,i1) - pasym(i2,i3,i1)
        end do
      end do
      
      do concurrent ( i2=1:8, i1=1:nback )
        psumN(i1,i2,2,0) = zero
        psumS(i1,i2,2,0) = zero
      end do
      
      !**************************************************************************************************************!
      !The backward (towards grid) fft, grid operations and the forward fft (towards space) *************************!
      !**************************************************************************************************************!
      call this%fourtrans%exec_c2r_sub( 8*nback, sumN(1) )
      call grid_sub( this%nFourier, sumN(1) )
      call this%fourtrans%exec_r2c_sub( 8*nforw, sumN(1) )
      
      call this%fourtrans%exec_c2r_sub( 8*nback, sumS(1) )
      call grid_sub( this%nFourier, sumS(1) )
      call this%fourtrans%exec_r2c_sub( 8*nforw, sumS(1) )
      
      !**************************************************************************************************************!
      !The forward (towards space) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      pssym(1:8,1:2,1:nforw) => ssym(1:16*nforw)
      pasym(1:8,1:2,1:nforw) => asym(1:16*nforw)
      
      psumN(1:nforw,1:8,1:2,0:this%jmax2) => sumN(1:16*nforw*this%jmax3)
      psumS(1:nforw,1:8,1:2,0:this%jmax2) => sumS(1:16*nforw*this%jmax3)
      
      psumN(1:nforw,1:8,2,0) = zero
      psumS(1:nforw,1:8,2,0) = zero
      
      do m = 0, this%jmax2
        do concurrent ( i1=1:nforw, i3=1:2, i2=1:8 )
          pssym(i2,i3,i1) = weight(i2) * ( psumN(i1,i2,i3,m) + psumS(i1,i2,i3,m) )
          pasym(i2,i3,i1) = weight(i2) * ( psumN(i1,i2,i3,m) - psumS(i1,i2,i3,m) )
        end do
        
        if ( m == 0 ) then
          pmm(1:8) = this%cmm(0)
        else
          pmm(1:8) = this%cmm(m) * sinx(1:8) * pmm(1:8)
        end if
        
        j = m
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          pmj2(1:8) = zero
          pmj1(1:8) = zero
          pmj(1:8)  = pmm(1:8)
          
          do concurrent ( i1=1:nforw, i3=1:2, i2=1:8 )
            pcr(i3,i1,mj) = pcr(i3,i1,mj) + pmj(i2) * pssym(i2,i3,i1)
          end do
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amj(mj-1) * cosx(i2) * pmj1(i2) - this%bmj(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw, i3=1:2, i2=1:8 )
            pcr(i3,i1,mj-1) = pcr(i3,i1,mj-1) + pmj(i2) * pasym(i2,i3,i1)
          end do
          
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amj(mj) * cosx(i2) * pmj1(i2) - this%bmj(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw, i3=1:2, i2=1:8 )
            pcr(i3,i1,mj) = pcr(i3,i1,mj) + pmj(i2) * pssym(i2,i3,i1)
          end do
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          pmj2(1:8) = pmj1(1:8)
          pmj1(1:8) = pmj(1:8)
          
          do concurrent ( i2=1:8 )
            pmj(i2)  = this%amj(mj+1) * cosx(i2) * pmj1(i2) - this%bmj(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw, i3=1:2, i2=1:8 )
            pcr(i3,i1,mj+1) = pcr(i3,i1,mj+1) + pmj(i2) * pasym(i2,i3,i1)
          end do
        end if
      end do
    end do
    
    !Cleaning
    deallocate( pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym, sumN, sumS )
    
  end subroutine lege_transform_sub
  
end submodule lege_transform