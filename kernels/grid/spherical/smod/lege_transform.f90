submodule (SphericalHarmonics) lege_transform
  implicit none; contains
  
  module subroutine lege_transform_sub(this, nforw, nback, cc, cr, grid_sub)
    class(T_lateralGrid),      intent(in)    :: this
    integer,                   intent(in)    :: nforw, nback
    complex(kind=dbl), target, intent(in)    :: cc(nback,*)
    complex(kind=dbl), target, intent(inout) :: cr(nforw,*)
    integer                                  :: i1, i2, i3, i, j, m, mj
    real(kind=dbl),            allocatable   :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), weight(:), grid(:)
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
    
    !Allocating needed memory :: step is set to 4
    allocate( pmm(4), pmj(4), pmj1(4), pmj2(4), cosx(4), sinx(4), weight(4), ssym(8*nback), asym(8*nback), &
            & sumN(8*nback*this%jmax3), sumS(8*nback*this%jmax3), grid(4*nback*this%nFourier) )
    
    !Stepping of the algorithm :: 4
    do i = 1, this%nLegendre, 4
      cosx(1:4)   = this%cosx(i:i+3)
      sinx(1:4)   = sqrt(1-cosx(1:4)**2)
      weight(1:4) = this%weight(i:i+3)
      
      sumN = zero
      sumS = zero
      
      !**************************************************************************************************************!
      !The backward (towards grid) sum over associated Legendre polynomials *****************************************!
      !**************************************************************************************************************!
      pssym(1:2,1:4,1:nback) => ssym(1:8*nback)
      pasym(1:2,1:4,1:nback) => asym(1:8*nback)
      
      psumN(1:2,1:nback,1:4,0:this%jmax2) => sumN(1:8*nback*this%jmax3)
      psumS(1:2,1:nback,1:4,0:this%jmax2) => sumS(1:8*nback*this%jmax3)
      
      do m = 0, this%jmax2
        ssym = zero
        asym = zero
        
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
          
          do concurrent ( i1 = 1:nback, i2 = 1:4, i3 = 1:2 )
            pssym(i3,i2,i1) = pssym(i3,i2,i1) + pcc(i3,i1,mj) * pmj(i2)
          end do
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amj(mj-1) * cosx(i2) * pmj1(i2) - this%bmj(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:4, i3=1:2 )
            pasym(i3,i2,i1) = pasym(i3,i2,i1) + pcc(i3,i1,mj-1) * pmj(i2)
          end do
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amj(mj) * cosx(i2) * pmj1(i2) - this%bmj(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:4, i3=1:2 )
            pssym(i3,i2,i1) = pssym(i3,i2,i1) + pcc(i3,i1,mj) * pmj(i2)
          end do
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amj(mj+1) * cosx(i2) * pmj1(i2) - this%bmj(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nback, i2=1:4, i3=1:2 )
            pasym(i3,i2,i1) = pasym(i3,i2,i1) + pcc(i3,i1,mj+1) * pmj(i2)
          end do
        end if
        
        do concurrent ( i2=1:4, i1=1:nback, i3=1:2 )
          psumN(i3,i1,i2,m) = pssym(i3,i2,i1) + pasym(i3,i2,i1)
          psumS(i3,i1,i2,m) = pssym(i3,i2,i1) - pasym(i3,i2,i1)
        end do
      end do
      
      do concurrent ( i2=1:4, i1=1:nback )
        psumN(2,i1,i2,0) = zero
        psumS(2,i1,i2,0) = zero
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
      pssym(1:2,1:4,1:nforw) => ssym(1:8*nforw)
      pasym(1:2,1:4,1:nforw) => asym(1:8*nforw)
      
      psumN(1:2,1:nforw,1:4,0:this%jmax2) => sumN(1:8*nforw*this%jmax3)
      psumS(1:2,1:nforw,1:4,0:this%jmax2) => sumS(1:8*nforw*this%jmax3)
      
      do m = 0, this%jmax2
        do concurrent ( i2=1:4, i1=1:nforw, i3=1:2 )
          pssym(i3,i2,i1) = weight(i2) * ( psumN(i3,i1,i2,m) + psumS(i3,i1,i2,m) )
          pasym(i3,i2,i1) = weight(i2) * ( psumN(i3,i1,i2,m) - psumS(i3,i1,i2,m) )
        end do
        
        if ( m == 0 ) then
          do concurrent ( i1=1:nforw, i2=1:4 )
            pssym(2,i2,i1) = zero
            pasym(2,i2,i1) = zero
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
          
          do concurrent ( i1=1:nforw , i2=1:4, i3=1:2 )
            pcr(i3,i1,mj) = pcr(i3,i1,mj) + pmj(i2) * pssym(i3,i2,i1)
          end do
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2)  = this%amj(mj-1) * cosx(i2) * pmj1(i2) - this%bmj(mj-1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:4, i3=1:2 )
            pcr(i3,i1,mj-1) = pcr(i3,i1,mj-1) + pmj(i2) * pasym(i3,i2,i1)
          end do
          
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2) = this%amj(mj) * cosx(i2) * pmj1(i2) - this%bmj(mj) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:4, i3=1:2 )
            pcr(i3,i1,mj) = pcr(i3,i1,mj) + pmj(i2) * pssym(i3,i2,i1)
          end do
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          pmj2(1:4) = pmj1(1:4)
          pmj1(1:4) = pmj(1:4)
          
          do concurrent ( i2=1:4 )
            pmj(i2) = this%amj(mj+1) * cosx(i2) * pmj1(i2) - this%bmj(mj+1) * pmj2(i2)
          end do
          
          do concurrent ( i1=1:nforw , i2=1:4, i3=1:2 )
            pcr(i3,i1,mj+1) = pcr(i3,i1,mj+1) + pmj(i2) * pasym(i3,i2,i1)
          end do
        end if
      end do
    end do
    
    !Cleaning
    deallocate( pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym, sumN, sumS, grid )
    
  end subroutine lege_transform_sub
  
end submodule lege_transform