submodule (SphericalHarmonics) vcvv
  implicit none ; contains
  
  module pure subroutine vcvv_sub(this, cajml, cbjml, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajml(:), cbjml(:)
    complex(kind=dbl),    intent(out) :: cjm(:)
    integer                           :: i, k, j, m, l, mj, lmj, i1, i2
    real(kind=dbl)                    :: cleb1, cleb2, cleb3, fac
    real(kind=dbl),     allocatable   :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), fftLege(:), fft(:,:), grid(:,:,:), sinx(:)
    complex(kind=dbl),  allocatable   :: cc(:,:), sumLegendreN(:,:,:), sumLegendreS(:,:,:), cr(:)
    complex(kind=dbl),  allocatable   :: symL(:,:), asymL(:,:), symF(:), asymF(:), sumFourierN(:,:), sumFourierS(:,:)
    
    allocate( cc(6,this%jms2) ); cc = czero
    
    mj = 0
      do m = 0, this%maxj
        do j = m, this%maxj
          mj = mj + 1
          
          if ( m == 0 ) then
            do l = abs(j-1), min(this%jmax, j+1)
              lmj = 3*(l*(l+1)/2+m+1)+j-l
              
              cleb1 = cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(l+j)
              cleb2 = cleb1_fn(j,m,1,+1,l,m+1)
              cleb3 = cleb1_fn(j,m,1, 0,l,m  )
              
              cc(1,mj) = cc(1,mj) + conjg( cajml(lmj  ) ) * cleb1
              cc(2,mj) = cc(2,mj) +        cajml(lmj  )   * cleb2
              cc(3,mj) = cc(3,mj) +        cajml(lmj-3)   * cleb3
              cc(4,mj) = cc(4,mj) + conjg( cbjml(lmj  ) ) * cleb1
              cc(5,mj) = cc(5,mj) +        cbjml(lmj  )   * cleb2
              cc(6,mj) = cc(6,mj) +        cbjml(lmj-3)   * cleb3
            end do
          else
            do l = abs(j-1), min(this%jmax, j+1)
              lmj = 3*(l*(l+1)/2+m-1)+j-l
              
              if (l > m+0) then
                cleb1 = cleb1_fn(j,m,1,-1,l,m-1)
                cleb2 = cleb1_fn(j,m,1,+1,l,m+1)
                cleb3 = cleb1_fn(j,m,1, 0,l,m  )
                
                cc(1,mj) = cc(1,mj) + cajml(lmj  ) * cleb1
                cc(2,mj) = cc(2,mj) + cajml(lmj+6) * cleb2
                cc(3,mj) = cc(3,mj) + cajml(lmj+3) * cleb3
                cc(4,mj) = cc(4,mj) + cbjml(lmj  ) * cleb1
                cc(5,mj) = cc(5,mj) + cbjml(lmj+6) * cleb2
                cc(6,mj) = cc(6,mj) + cbjml(lmj+3) * cleb3
                
              else if (l > m-1) then
                cleb1 = cleb1_fn(j,m,1,-1,l,m-1)
                cleb3 = cleb1_fn(j,m,1, 0,l,m  )
                
                cc(1,mj) = cc(1,mj) + cajml(lmj  ) * cleb1
                cc(3,mj) = cc(3,mj) + cajml(lmj+3) * cleb3
                cc(4,mj) = cc(4,mj) + cbjml(lmj  ) * cleb1
                cc(6,mj) = cc(6,mj) + cbjml(lmj+3) * cleb3
                
              else
                cleb1 = cleb1_fn(j,m,1,-1,l,m-1)
                
                cc(1,mj) = cc(1,mj) + cajml(lmj) * cleb1
                cc(4,mj) = cc(4,mj) + cbjml(lmj) * cleb1
              end if
            end do
          end if
          
          cc(1,mj) =           cc(1,mj)-cc(2,mj)          
          cc(2,mj) = cunit * ( cc(1,mj)+cc(2,mj)+cc(2,mj) )
          cc(3,mj) =     2 *   cc(3,mj)               
          cc(4,mj) =           cc(4,mj)-cc(5,mj)          
          cc(5,mj) = cunit * ( cc(4,mj)+cc(5,mj)+cc(5,mj) )
        end do
      end do
    
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), sinx(step), fftLege(step), symL(step,6), asymL(step,6),  & 
            & sumLegendreN(6,step,0:this%maxj), sumLegendreS(6,step,0:this%maxj), grid(6,step,0:this%nFourier-1), cr(this%jms1), &
            & fft(step,0:this%nFourier-1), sumFourierN(step,0:this%jmax+1), sumFourierS(step,0:this%jmax+1), symF(step),         &
            & asymF(step) ) ; cr = czero
    
    do i = 1, this%nLegendre, step
      cosx    = this%roots(i:i+step-1)
      sinx    = sqrt(1-cosx**2)
      fftLege = this%fftLege(i:i+step-1)
      
      pmm = 1._dbl
      
      sumLegendreN = czero
      sumLegendreS = czero
      
      m = 0
        pmj2 = 0._dbl
        pmj1 = 0._dbl
        pmj  = pmm
        
        symL  = czero
        asymL = czero
        
        j = m
          mj = m*(this%maxj+1)-m*(m+1)/2+j+1
          
          do concurrent ( i1=1:6 , i2=1:step )
            symL(i2,i1) = symL(i2,i1) + cc(i1,mj)
          end do
          
        do j = 1, (this%maxj-m)/2
          mj = mj+2
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj-1) * cosx * pmj1 - this%bmjrr(mj-1) * pmj2
          
          do concurrent ( i1=1:6 , i2=1:step )
            asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj-1) * pmj(i2)
          end do
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
          
          do concurrent ( i1=1:6 , i2=1:step )
            symL(i2,i1) = symL(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
          
          if ( maxval(abs(pmj)) < this%tolm ) exit
        end do
        
        if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
          mj = mj+1
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
          
          do concurrent ( i1=1:6 , i2=1:step )
            asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        end if
        
        do concurrent ( i2=1:step, i1=1:6 )
          sumLegendreN(i1,i2,0)%re = symL(i2,i1)%re + asymL(i2,i1)%re ; sumLegendreN(i1,i2,0)%im = 0._dbl
          sumLegendreS(i1,i2,0)%re = symL(i2,i1)%re - asymL(i2,i1)%re ; sumLegendreS(i1,i2,0)%im = 0._dbl
        end do
        
      do m = 1, this%maxj
        fac = -sqrt( ( 2*m+1 ) / ( 2._dbl*m ) ) ; pmm = fac * sinx * pmm ; if (maxval(abs(pmm)) < this%tolm) exit
        
        pmj2 = 0._dbl
        pmj1 = 0._dbl
        pmj  = pmm
        
        symL  = czero
        asymL = czero
        
        j = m
          mj = m*(this%maxj+1)-m*(m+1)/2+j+1
          
          do concurrent ( i1=1:6 , i2=1:step )
            symL(i2,i1) = symL(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        
        do j = 1, (this%maxj-m)/2
          mj = mj+2
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj-1) * cosx * pmj1 - this%bmjrr(mj-1) * pmj2
          
          do concurrent ( i1=1:6 , i2=1:step )
            asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj-1) * pmj(i2)
          end do
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj  = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
          
          do concurrent ( i1=1:6 , i2=1:step )
            symL(i2,i1) = symL(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
          
          if ( maxval(abs(pmj)) < this%tolm ) exit
        end do
        
        if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
          mj = mj+1
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj  = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
          
          do concurrent ( i1=1:6 , i2=1:step )
            asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        end if
        
        do concurrent ( i2=1:step, i1=1:6 )
          sumLegendreN(i1,i2,m) = ( symL(i2,i1) + asymL(i2,i1) )
          sumLegendreS(i1,i2,m) = ( symL(i2,i1) - asymL(i2,i1) ) 
        end do
      end do
      
      call this%fourtrans%exec_c2r_sub(6*step, this%maxj+1, sumLegendreN, grid)
      
      do concurrent ( i1=0:this%nFourier-1, i2=1:step )
        fft(i2,i1) = sum( grid(1:3,i2,i1) * grid(4:6,i2,i1) )
      end do
      
      call this%fourtrans%exec_r2c_sub(step, this%maxj, fft, sumFourierN)
      
      call this%fourtrans%exec_c2r_sub(6*step, this%maxj+1, sumLegendreS, grid)
      
      do concurrent ( i1=0:this%nFourier-1, i2=1:step )
        fft(i2,i1) = sum( grid(1:3,i2,i1) * grid(4:6,i2,i1) )
      end do
      
      call this%fourtrans%exec_r2c_sub(step, this%maxj, fft, sumFourierS)
      
      pmm  = 1._dbl
      
      m = 0
        pmj2 = 0._dbl
        pmj1 = 0._dbl
        pmj  = pmm
        
        do concurrent ( i2=1:step )
          symF(i2)  = fftLege(i2) * ( sumFourierN(i2,m) + sumFourierS(i2,m) ) ; symF(i2)%im  = 0._dbl
          asymF(i2) = fftLege(i2) * ( sumFourierN(i2,m) - sumFourierS(i2,m) ) ; asymF(i2)%im = 0._dbl
        end do
        
        j = m
          mj = m*this%maxj-m*(m+1)/2+j+1
          
          cr(mj) = cr(mj) + sum( symF(:) )
            
        do j = 1, (this%jmax-m+1)/2
          mj = mj+2
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj-1+m) * cosx * pmj1 - this%bmjrr(mj-1+m) * pmj2
          
          cr(mj-1) = cr(mj-1) + sum( pmj(:) * asymF(:) )
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
          
          cr(mj) = cr(mj) + sum( pmj(:) * symF(:) )
        end do
        
        if (mod(this%jmax+1-m,2) /= 0) then
          mj = mj+1
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
          
          cr(mj) = cr(mj) + sum( pmj(:) * asymF(:) )          
        end if
        
      do m = 1, this%jmax+1
        fac = -sqrt( ( 2*m+1 ) / ( 2._dbl*m ) ) ; pmm = fac * sinx * pmm
        
        pmj2 = 0._dbl
        pmj1 = 0._dbl
        pmj  = pmm
        
        do concurrent ( i2=1:step )
          symF(i2)  = fftLege(i2) * ( sumFourierN(i2,m) + sumFourierS(i2,m) )
          asymF(i2) = fftLege(i2) * ( sumFourierN(i2,m) - sumFourierS(i2,m) )
        end do
        
        j = m
          mj = m*this%maxj-m*(m+1)/2+j+1
          
          cr(mj) = cr(mj) + sum( pmj(:) * symF(:) )
          
        do j = 1, (this%jmax+1-m)/2
          mj = mj+2
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj-1+m) * cosx * pmj1 - this%bmjrr(mj-1+m) * pmj2
          
          cr(mj-1) = cr(mj-1) + sum( pmj(:) * asymF(:) )
                    
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
          
          cr(mj) = cr(mj) + sum( pmj(:) * symF(:) )
        end do
        
        if (mod(this%jmax+1-m,2) /= 0) then
          mj = mj+1
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
          
          cr(mj) = cr(mj) + sum( pmj(:) * asymF(:) )
        end if
      end do
    end do
    
    deallocate( cc, fft, sumFourierN, sumFourierS, pmm, pmj, pmj1, pmj2, cosx, sinx, sumLegendreN, sumLegendreS, &
              & fftLege, grid, symL, asymL, symF, asymF )
    
    fac = 1 / (16 * this%nLegendre**2 * sqrt(pi) )
    
    do j = 0, this%jmax
      m = 0
        cjm(j*(j+1)/2+m+1)%re = cr(m*(this%maxj)-m*(m+1)/2+j+1)%re * fac
        cjm(j*(j+1)/2+m+1)%im = 0._dbl
      
      do m = 1, j
        cjm(j*(j+1)/2+m+1) = cr(m*(this%maxj)-m*(m+1)/2+j+1) * fac
      end do
    end do
    
    deallocate(cr)
    
  end subroutine vcvv_sub
  
end submodule vcvv