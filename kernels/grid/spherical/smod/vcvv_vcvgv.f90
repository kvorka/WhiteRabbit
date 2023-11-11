submodule (SphericalHarmonics) vcvv_vcvgv
  implicit none
  
  contains
  
  module pure subroutine vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(:), q(:), v(:)
    complex(kind=dbl),    intent(out) :: cjm(:,:)
    integer                           :: i, j, m, l, ijm, ijml, lmj, mj, mj1, mj2, i1, i2
    real(kind=dbl)                    :: cleb, fac1, fac2
    real(kind=dbl),       allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), fftLege(:), grid(:,:,:), fft(:,:,:)
    complex(kind=dbl)                 :: cr12
    complex(kind=dbl),    allocatable :: sum1(:), sum2(:), sum3(:)
    complex(kind=dbl),    allocatable :: cab(:,:), cc(:,:), cr(:,:), symL(:,:), asymL(:,:), symF(:,:), asymF(:,:)
    complex(kind=dbl),    allocatable :: sumLegendreN(:,:,:), sumLegendreS(:,:,:), sumFourierN(:,:,:), sumFourierS(:,:,:)
    
    allocate( cab(5,this%jmv1), sum1(2), sum2(2), sum3(2) ) ; cab = czero
    
    do concurrent( ijml = 1:this%jmv )
      cab(1,ijml) = q(ijml)
      cab(2,ijml) = v(ijml)
    end do
    
    do j = 0, this%jmax+1
      fac1 = sqrt((j  )/(2*j+1._dbl))
      fac2 = sqrt((j+1)/(2*j+1._dbl))
      
      do m = 0, j
        sum1 = czero ; sum2 = czero ; sum3 = czero
        
        if (m == 0) then
          do l = abs(j-1), min(this%jmax, j+1)
            lmj = 3*(l*(l+1)/2+m+1)+j-l
            
            cleb = cleb1_fn(j,m,1,0,l,m)
              sum3(1) = sum3(1) + ( dv_r(lmj-3) + (j+1) * v(lmj-3) / ri ) * cleb
              sum3(2) = sum3(2) - ( dv_r(lmj-3) - (j  ) * v(lmj-3) / ri ) * cleb
            
            cleb = cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(j+l)
              sum1(1) = sum1(1) + conjg( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
              sum1(2) = sum1(2) - conjg( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
            
            cleb = cleb1_fn(j,m,1,+1,l,m+1)
              sum2(1) = sum2(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
              sum2(2) = sum2(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
          end do
        
        else
          do l = max(abs(m-1), j-1), min(this%jmax, j+1)
            lmj = 3*(l*(l+1)/2+m-1)+j-l
            
            if (l > m) then
              cleb = cleb1_fn(j,m,1,-1,l,m-1)
                sum1(1) = sum1(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
                sum1(2) = sum1(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
              
              cleb = cleb1_fn(j,m,1,0,l,m)
                sum3(1) = sum3(1) + ( dv_r(lmj+3) + (j+1) * v(lmj+3) / ri ) * cleb
                sum3(2) = sum3(2) - ( dv_r(lmj+3) - (j  ) * v(lmj+3) / ri ) * cleb
              
              cleb = cleb1_fn(j,m,1,+1,l,m+1)
                sum2(1) = sum2(1) + ( dv_r(lmj+6) + (j+1) * v(lmj+6) / ri ) * cleb
                sum2(2) = sum2(2) - ( dv_r(lmj+6) - (j  ) * v(lmj+6) / ri ) * cleb
            
            else if (l > m-1) then
              cleb = cleb1_fn(j,m,1,-1,l,m-1)
                sum1(1) = sum1(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
                sum1(2) = sum1(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
              
              cleb = cleb1_fn(j,m,1,0,l,m)
                sum3(1) = sum3(1) + ( dv_r(lmj+3) + (j+1) * v(lmj+3) / ri ) * cleb
                sum3(2) = sum3(2) - ( dv_r(lmj+3) - (j  ) * v(lmj+3) / ri ) * cleb
              
            else
              cleb = cleb1_fn(j,m,1,-1,l,m-1)
                sum1(1) = sum1(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
                sum1(2) = sum1(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
            end if
          end do
        end if
        
        ijml = 3*(j*(j+1)/2+m)+abs(j-1)-j
          cab(3,ijml) =         ( sum1(1) - sum2(1) ) * fac1
          cab(4,ijml) = cunit * ( sum1(1) + sum2(1) ) * fac1
          cab(5,ijml) =         ( sum3(1)           ) * fac1
        
        ijml = 3*(j*(j+1)/2+m)+(j+1)-j
          cab(3,ijml) =         ( sum1(2) - sum2(2) ) * fac2
          cab(4,ijml) = cunit * ( sum1(2) + sum2(2) ) * fac2
          cab(5,ijml) =         ( sum3(2)           ) * fac2
      end do
    end do
    
    deallocate( sum1, sum2, sum3 )
    allocate( cc(15, this%jms2), sum1(5), sum2(5), sum3(5) ) ; cc = czero
    
    mj = 0
      do m = 0, this%maxj
        do j = m, this%maxj
          sum1 = czero ; sum2 = czero ; sum3 = czero
          
          if (m == 0) then
            do l = abs(j-1), min(this%jmax+1, j+1)
              lmj = 3*(l*(l+1)/2+m+1)+j-l
              
              cleb = cleb1_fn(j,m,1,0,l,m)
                do concurrent ( i1 = 1:5 )
                  sum3(i1) = sum3(i1) + cab(i1,lmj-3) * cleb
                end do
              
              cleb = cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(j+l)
                do concurrent ( i1 = 1:5 )
                  sum1(i1) = sum1(i1) + conjg( cab(i1,lmj) ) * cleb
                end do
              
              cleb = cleb1_fn(j,m,1,+1,l,m+1)
                do concurrent ( i1 = 1:5 )
                  sum2(i1) = sum2(i1) + cab(i1,lmj) * cleb
                end do
            end do
          else
            do l = abs(j-1), min(this%jmax+1, j+1)
              lmj = 3*(l*(l+1)/2+m-1)+j-l
              
              if (l > m) then
                cleb = cleb1_fn(j,m,1,-1,l,m-1)
                  do concurrent ( i1 = 1:5 )
                    sum1(i1) = sum1(i1) + cab(i1,lmj) * cleb
                  end do
                
                cleb = cleb1_fn(j,m,1,0,l,m)
                  do concurrent ( i1 = 1:5 )
                    sum3(i1) = sum3(i1) + cab(i1,lmj+3) * cleb
                  end do
                
                cleb = cleb1_fn(j,m,1,+1,l,m+1)
                  do concurrent ( i1 = 1:5 )
                    sum2(i1) = sum2(i1) + cab(i1,lmj+6) * cleb
                  end do
                  
              else if (l > m-1) then
                cleb = cleb1_fn(j,m,1,-1,l,m-1)
                  do concurrent ( i1 = 1:5 )
                    sum1(i1) = sum1(i1) + cab(i1,lmj) * cleb
                  end do
                
                cleb = cleb1_fn(j,m,1,0,l,m)
                  do concurrent ( i1 = 1:5 )
                    sum3(i1) = sum3(i1) + cab(i1,lmj+3) * cleb
                  end do
              
              else
                cleb = cleb1_fn(j,m,1,-1,l,m-1)
                  do concurrent (i1 = 1:5)
                    sum1(i1) = sum1(i1) + cab(i1,lmj) * cleb
                  end do
              end if
            end do
          end if
          
          mj = mj+1
            do concurrent ( i1 = 1:5 )
              i2 = 3*(i1-1)+1
              
              cc(i2  ,mj) =           sum1(i1) - sum2(i1)
              cc(i2+1,mj) = cunit * ( sum1(i1) + sum2(i1) )
              cc(i2+2,mj) =           sum3(i1)
            end do
        end do
      end do
      
      do concurrent ( mj = 1:this%jms2 )
        cc(1: 5,mj) = cc(1: 5,mj) / 2
        cc(7:15,mj) = cc(7:15,mj) / 2
      end do
    
    deallocate(cab, sum1, sum2, sum3)
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), fftLege(step), sinx(step), symL(step,15), asymL(step,15), &
            & sumLegendreN(15,step,0:this%maxj), sumLegendreS(15,step,0:this%maxj), grid(15,step,0:this%nFourier-1),              &
            & cr(4,this%jms1), fft(4,step,0:this%nFourier-1), sumFourierN(4,step,0:this%jmax+1), symF(step,4), asymF(step,4),     &
            & sumFourierS(4,step,0:this%jmax+1) ) ; cr = czero
    
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
          
          do concurrent ( i1=1:15 , i2=1:step )
            symL(i2,i1) = symL(i2,i1) + cc(i1,mj)
          end do
          
        do j = 1, (this%maxj-m)/2
          mj = mj+2
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj-1) * cosx * pmj1 - this%bmjrr(mj-1) * pmj2
          
          do concurrent ( i1=1:15 , i2=1:step )
            asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj-1) * pmj(i2)
          end do
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
          
          do concurrent ( i1=1:15 , i2=1:step )
            symL(i2,i1) = symL(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
          
          if ( maxval(abs(pmj)) < this%tolm ) exit
        end do
        
        if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
          mj = mj+1
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
          
          do concurrent ( i1=1:15 , i2=1:step )
            asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        end if
        
        do concurrent ( i2=1:step, i1=1:15 )
          sumLegendreN(i1,i2,0)%re = symL(i2,i1)%re + asymL(i2,i1)%re ; sumLegendreN(i1,i2,0)%im = 0._dbl
          sumLegendreS(i1,i2,0)%re = symL(i2,i1)%re - asymL(i2,i1)%re ; sumLegendreS(i1,i2,0)%im = 0._dbl
        end do
        
      do m = 1, this%maxj
        fac1 = -sqrt( ( 2*m+1 ) / ( 2._dbl*m ) ) ; pmm = fac1 * sinx * pmm ; if (maxval(abs(pmm)) < this%tolm) exit
        
        pmj2 = 0._dbl
        pmj1 = 0._dbl
        pmj  = pmm
        
        symL  = czero
        asymL = czero
        
        j = m
          mj = m*(this%maxj+1)-m*(m+1)/2+j+1
          
          do concurrent ( i1=1:15 , i2=1:step )
            symL(i2,i1) = symL(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        
        do j = 1, (this%maxj-m)/2
          mj = mj+2
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj-1) * cosx * pmj1 - this%bmjrr(mj-1) * pmj2
          
          do concurrent ( i1=1:15 , i2=1:step )
            asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj-1) * pmj(i2)
          end do
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj  = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
          
          do concurrent ( i1=1:15 , i2=1:step )
            symL(i2,i1) = symL(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
          
          if ( maxval(abs(pmj)) < this%tolm ) exit
        end do
        
        if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
          mj = mj+1
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj  = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
          
          do concurrent ( i1=1:15 , i2=1:step )
            asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj) * pmj(i2)
          end do
        end if
        
        do concurrent ( i2=1:step, i1=1:15 )
          sumLegendreN(i1,i2,m) = ( symL(i2,i1) + asymL(i2,i1) )
          sumLegendreS(i1,i2,m) = ( symL(i2,i1) - asymL(i2,i1) ) 
        end do
      end do
      
      call this%fourtrans%exec_c2r_sub(15*step, this%maxj+1, sumLegendreN, grid)
      
      do concurrent (i1=0:this%nFourier-1, i2=1:step)
        fft(1,i2,i1) = grid( 1,i2,i1) * grid( 4,i2,i1) + &    
                     & grid( 2,i2,i1) * grid( 5,i2,i1) + &    
                     & grid( 3,i2,i1) * grid( 6,i2,i1)
        fft(2,i2,i1) = grid( 4,i2,i1) * grid( 7,i2,i1) + &
                     & grid( 5,i2,i1) * grid( 8,i2,i1) + &
                     & grid( 6,i2,i1) * grid( 9,i2,i1)
        fft(3,i2,i1) = grid( 4,i2,i1) * grid(10,i2,i1) + &
                     & grid( 5,i2,i1) * grid(11,i2,i1) + &
                     & grid( 6,i2,i1) * grid(12,i2,i1)
        fft(4,i2,i1) = grid( 4,i2,i1) * grid(13,i2,i1) + &
                     & grid( 5,i2,i1) * grid(14,i2,i1) + &
                     & grid( 6,i2,i1) * grid(15,i2,i1)
      end do
      
      call this%fourtrans%exec_r2c_sub(4*step, this%maxj, fft, sumFourierN)
      
      call this%fourtrans%exec_c2r_sub(15*step, this%maxj+1, sumLegendreS, grid)
      
      do concurrent (i1=0:this%nFourier-1, i2=1:step)
        fft(1,i2,i1) = grid( 1,i2,i1) * grid( 4,i2,i1) + &    
                     & grid( 2,i2,i1) * grid( 5,i2,i1) + &    
                     & grid( 3,i2,i1) * grid( 6,i2,i1)
        fft(2,i2,i1) = grid( 4,i2,i1) * grid( 7,i2,i1) + &
                     & grid( 5,i2,i1) * grid( 8,i2,i1) + &
                     & grid( 6,i2,i1) * grid( 9,i2,i1)
        fft(3,i2,i1) = grid( 4,i2,i1) * grid(10,i2,i1) + &
                     & grid( 5,i2,i1) * grid(11,i2,i1) + &
                     & grid( 6,i2,i1) * grid(12,i2,i1)
        fft(4,i2,i1) = grid( 4,i2,i1) * grid(13,i2,i1) + &
                     & grid( 5,i2,i1) * grid(14,i2,i1) + &
                     & grid( 6,i2,i1) * grid(15,i2,i1)
      end do
      
      call this%fourtrans%exec_r2c_sub(4*step, this%maxj, fft, sumFourierS)
      
      pmm  = 1._dbl
      
      m = 0
        pmj2 = 0._dbl
        pmj1 = 0._dbl
        pmj  = pmm
        
        do concurrent (i2=1:step, i1=1:4)
          symF(i2,i1)  = fftLege(i2) * ( sumFourierN(i1,i2,m) + sumFourierS(i1,i2,m) ) ; symF(i2,i1)%im  = 0._dbl
          asymF(i2,i1) = fftLege(i2) * ( sumFourierN(i1,i2,m) - sumFourierS(i1,i2,m) ) ; asymF(i2,i1)%im = 0._dbl
        end do
        
        j = m
          mj = m*this%maxj-m*(m+1)/2+j+1
          
          do concurrent ( i1=1:4 , i2=1:step )
            cr(i1,mj) = cr(i1,mj) + symF(i2,i1)
          end do
        
        do j = 1, (this%jmax-m+1)/2
          mj = mj+2
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj-1+m) * cosx * pmj1 - this%bmjrr(mj-1+m) * pmj2
          
          do concurrent ( i1=1:4 , i2=1:step )
            cr(i1,mj-1) = cr(i1,mj-1) + pmj(i2) * asymF(i2,i1)
          end do
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
          
          do concurrent ( i1=1:4 , i2=1:step )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * symF(i2,i1)
          end do
        end do
        
        if (mod(this%jmax+1-m,2) /= 0) then
          mj = mj+1
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
          
          do concurrent ( i1=1:4 , i2=1:step )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * asymF(i2,i1)
          end do
        end if
      
      do m = 1, this%jmax+1
        fac1 = -sqrt( ( 2*m+1 ) / ( 2._dbl*m ) ) ; pmm = fac1 * sinx * pmm
        
        pmj2 = 0._dbl
        pmj1 = 0._dbl
        pmj  = pmm
        
        do concurrent ( i2=1:step, i1=1:4 )
          symF(i2,i1)  = fftLege(i2) * ( sumFourierN(i1,i2,m) + sumFourierS(i1,i2,m) )
          asymF(i2,i1) = fftLege(i2) * ( sumFourierN(i1,i2,m) - sumFourierS(i1,i2,m) )
        end do
        
        j = m
          mj = m*this%maxj-m*(m+1)/2+j+1
          
          do concurrent ( i1=1:4 , i2=1:step )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * symF(i2,i1)
          end do
        
        do j = 1, (this%jmax+1-m)/2
          mj = mj+2
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj-1+m) * cosx * pmj1 - this%bmjrr(mj-1+m) * pmj2
          
          do concurrent ( i1=1:4 , i2=1:step )
            cr(i1,mj-1) = cr(i1,mj-1) + pmj(i2) * asymF(i2,i1)
          end do
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
          
          do concurrent ( i1=1:4 , i2=1:step )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * symF(i2,i1)
          end do
        end do
        
        if (mod(this%jmax+1-m,2) /= 0) then
          mj = mj+1
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
          
          do concurrent ( i1=1:4 , i2=1:step )
            cr(i1,mj) = cr(i1,mj) + pmj(i2) * asymF(i2,i1)
          end do
        end if
      end do
    end do
    
    deallocate( cc, sumLegendreN, sumLegendreS, grid, fft, sumFourierN, sumFourierS, pmm, pmj, pmj1, pmj2, cosx, sinx, fftLege, &
              & symL, asymL, symF, asymF )
    
    fac1 = 1 / ( 4 * this%nLegendre**2 * sqrt(pi) )
    
    do concurrent ( mj = 1:this%jms1 ) 
      cr(1,mj) = cr(1,mj) * fac1
      cr(2,mj) = cr(2,mj) * fac1 / 2
      cr(3,mj) = cr(3,mj) * fac1 / 2 * cunit
      cr(4,mj) = cr(4,mj) * fac1
      
      cr12     = cr(2,mj) - cr(3,mj)
      cr(3,mj) = cr(2,mj) + cr(3,mj)
      cr(2,mj) = cr12
    end do
    
    j = 0
      m = 0
        ijm = 1
        mj  = m*(this%maxj)-m*(m+1)/2+j+1
        mj2 = mj + this%maxj - m
        
        cjm(1,ijm) =        cr(1,mj  )                                ; cjm(1,ijm)%im = 0._dbl
        cjm(4,ijm) =        cr(2,mj2 )   * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                   &        cr(4,mj+1)   * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                   & conjg( cr(2,mj2 ) ) * cleb1_fn(j+1,m-1,1,+1,j,m) ; cjm(4,ijm)%im = 0._dbl
    
    do j = 1, this%jmax
      m = 0
        ijm = ijm+1
        mj  = m*(this%maxj)-m*(m+1)/2+j+1
        mj2 = mj + this%maxj - m - 1
        
        cjm(1,ijm) =        cr(1,mj   )                                ; cjm(1,ijm)%im = 0._dbl
        cjm(2,ijm) =        cr(2,mj2-1)   * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                   &        cr(4,mj -1)   * cleb1_fn(j-1,m+0,1, 0,j,m) + &
                   & conjg( cr(2,mj2-1) ) * cleb1_fn(j-1,m-1,1,+1,j,m) ; cjm(2,ijm)%im = 0._dbl
        cjm(3,ijm) =        cr(2,mj2  )   * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                   &        cr(4,mj   )   * cleb1_fn(j  ,m+0,1, 0,j,m) + &
                   & conjg( cr(2,mj2  ) ) * cleb1_fn(j  ,m-1,1,+1,j,m) ; cjm(3,ijm)%re = 0._dbl
        cjm(4,ijm) =        cr(2,mj2+1)   * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                   &        cr(4,mj +1)   * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                   & conjg( cr(2,mj2+1) ) * cleb1_fn(j+1,m-1,1,+1,j,m) ; cjm(4,ijm)%im = 0._dbl
      
      do m = 1, j
        ijm = ijm+1
        mj  = m*(this%maxj)-m*(m+1)/2+j+1
        mj1 = mj - this%maxj + m
        mj2 = mj + this%maxj - m - 1
        
        cjm(1,ijm) = cr(1,mj   )
        cjm(2,ijm) = cr(2,mj2-1) * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                   & cr(4,mj -1) * cleb1_fn(j-1,m+0,1, 0,j,m) - &
                   & cr(3,mj1-1) * cleb1_fn(j-1,m-1,1,+1,j,m)
        cjm(3,ijm) = cr(2,mj2  ) * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                   & cr(4,mj   ) * cleb1_fn(j  ,m+0,1, 0,j,m) - &
                   & cr(3,mj1  ) * cleb1_fn(j  ,m-1,1,+1,j,m)
        cjm(4,ijm) = cr(2,mj2+1) * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                   & cr(4,mj +1) * cleb1_fn(j+1,m+0,1, 0,j,m) - &
                   & cr(3,mj1+1) * cleb1_fn(j+1,m-1,1,+1,j,m)
      end do
    end do
    
    deallocate(cr)
    
  end subroutine vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv