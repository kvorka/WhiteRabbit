submodule (SphericalHarmonics) vcsv_vcvv_vcvgv
  implicit none

  contains

  subroutine init_vcsv_vcvv_vcvgv_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    
    call this%fourtrans%init_sub( this%nFourier )
    
  end subroutine init_vcsv_vcvv_vcvgv_sub

  pure subroutine vcsv_vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(:), q(:), v(:)
    complex(kind=dbl),    intent(out) :: cjm(:,:)
    integer                           :: i, j, m, l, ijm, lm, ijml, i1, i2, mj, s, iL, lm1, lm2
    real(kind=dbl)                    :: fac, cleb1, cleb2, cleb3
    real(kind=dbl),       allocatable :: gc(:), pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), fftLege(:)
    real(kind=dbl),       allocatable :: grid(:,:,:), fft(:,:,:)
    complex(kind=dbl)                 :: cr12
    complex(kind=dbl),    allocatable :: sum1(:), sum2(:), sum3(:)
    complex(kind=dbl),    allocatable :: cab(:,:), cc(:,:), cr(:,:), symL(:,:), asymL(:,:), symF(:,:), asymF(:,:)
    complex(kind=dbl),    allocatable :: sumLegendreN(:,:,:), sumLegendreS(:,:,:), sumFourierN(:,:,:), sumFourierS(:,:,:)

    allocate(cab(6,this%jmv1), sum1(2), sum2(2), sum3(2), gc(2)) ; cab = czero

    do ijml = 1, this%jmv
      cab(1,ijml) = q(ijml)
      cab(2,ijml) = v(ijml)
      cab(6,ijml) = dv_r(ijml)
    end do

    do j = 1, this%jmax+1
      gc(1) = sqrt(j*(j+1)*(j+1)/(2*j+1._dbl))
      gc(2) = sqrt(j*(j  )*(j+1)/(2*j+1._dbl))

      do m = 0, j
        sum1 = czero
        sum2 = czero
        sum3 = czero

        if (m == 0) then
          do l = abs(j-1), min(this%jmax, j+1)
            lm = 3*(l*(l+1)/2+m+1)+j-l

            cleb1 = cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(j+l)
            cleb3 = cleb1_fn(j,m,1, 0,l,m+0)
            cleb2 = cleb1_fn(j,m,1,+1,l,m+1)

            do concurrent ( i1=1:2 )
              sum1(i1) = sum1(i1) + conjg( v(lm  ) ) * cleb1 * gc(i1)
              sum3(i1) = sum3(i1) +        v(lm-3)   * cleb3 * gc(i1)
              sum2(i1) = sum2(i1) +        v(lm  )   * cleb2 * gc(i1)
            end do
          end do
        else
          do l = max(abs(m-1), j-1), min(this%jmax, j+1)
            lm = 3*(l*(l+1)/2+m-1)+j-l

            if (l > m) then
              cleb1 = cleb1_fn(j,m,1,-1,l,m-1)
              cleb3 = cleb1_fn(j,m,1, 0,l,m+0)
              cleb2 = cleb1_fn(j,m,1,+1,l,m+1)

              do concurrent ( i1=1:2 )
                sum1(i1) = sum1(i1) + v(lm  ) * cleb1 * gc(i1)
                sum3(i1) = sum3(i1) + v(lm+3) * cleb3 * gc(i1)
                sum2(i1) = sum2(i1) + v(lm+6) * cleb2 * gc(i1)
              end do
            
            else if (l > m-1) then
              cleb1 = cleb1_fn(j,m,1,-1,l,m-1)
              cleb3 = cleb1_fn(j,m,1, 0,l,m+0)

              do concurrent ( i1=1:2 )
                sum1(i1) = sum1(i1) + v(lm  ) * cleb1 * gc(i1)
                sum3(i1) = sum3(i1) + v(lm+3) * cleb3 * gc(i1)
              end do
            
            else
              cleb1 = cleb1_fn(j,m,1,-1,l,m-1)

              do concurrent ( i1=1:2 )
                sum1(i1) = sum1(i1) + v(lm  ) * cleb1 * gc(i1)
              end do
            end if
          end do
        end if

        ijml = 3*(j*(j+1)/2+m)+abs(j-1)-j
          cab(3,ijml) =         ( sum1(1) - sum2(1) ) / ri
          cab(4,ijml) = cunit * ( sum1(1) + sum2(1) ) / ri
          cab(5,ijml) =         ( sum3(1)           ) / ri

        ijml = 3*(j*(j+1)/2+m)+(j+1)-j
          cab(3,ijml) =         ( sum1(2) - sum2(2) ) / ri
          cab(4,ijml) = cunit * ( sum1(2) + sum2(2) ) / ri
          cab(5,ijml) =         ( sum3(2)           ) / ri
      end do
    end do

    deallocate(sum1, sum2, sum3, gc)
    allocate(cc(19, this%jms2), sum1(6), sum2(6), sum3(6)) ; cc = czero
  
    mj = 0
      do m = 0, this%maxj
        do j = m, this%maxj
          sum1 = czero
          sum2 = czero
          sum3 = czero
    
          if (m == 0) then
            do l = abs(j-1), min(this%jmax+1, j+1)
              lm = 3*(l*(l+1)/2+m+1)+j-l
    
              cleb1 = cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(j+l)
              cleb3 = cleb1_fn(j,m,1, 0,l,m  )
              cleb2 = cleb1_fn(j,m,1,+1,l,m+1)
    
              do concurrent ( i1=1:6 )
                sum1(i1) = sum1(i1) + conjg( cab(i1,lm  ) ) * cleb1
                sum3(i1) = sum3(i1) +        cab(i1,lm-3)   * cleb3
                sum2(i1) = sum2(i1) +        cab(i1,lm  )   * cleb2
              end do
            end do
          else
            do l = abs(j-1), min(this%jmax+1, j+1)
              lm = 3*(l*(l+1)/2+m-1)+j-l
    
              if (l > m) then
                cleb1 = cleb1_fn(j,m,1,-1,l,m-1)
                cleb3 = cleb1_fn(j,m,1, 0,l,m  )
                cleb2 = cleb1_fn(j,m,1,+1,l,m+1)
    
                do concurrent (i1 = 1:6)
                  sum1(i1) = sum1(i1) + cab(i1,lm  ) * cleb1
                  sum3(i1) = sum3(i1) + cab(i1,lm+3) * cleb3
                  sum2(i1) = sum2(i1) + cab(i1,lm+6) * cleb2
                end do
    
              else if (l > m-1) then
                cleb1 = cleb1_fn(j,m,1,-1,l,m-1)
                cleb3 = cleb1_fn(j,m,1, 0,l,m  )
    
                do concurrent (i1 = 1:6)
                  sum1(i1) = sum1(i1) + cab(i1,lm  ) * cleb1
                  sum3(i1) = sum3(i1) + cab(i1,lm+3) * cleb3
                end do
              else
                cleb1 = cleb1_fn(j,m,1,-1,l,m-1)
    
                do concurrent (i1 = 1:6)
                  sum1(i1) = sum1(i1) + cab(i1,lm) * cleb1
                end do
              end if
            end do
          end if
    
          mj = mj+1
            do i1 = 1, 6
              i2 = 3*(i1-1)+1
    
              cc(i2  ,mj) =           sum1(i1) - sum2(i1)
              cc(i2+1,mj) = cunit * ( sum1(i1) + sum2(i1) )
              cc(i2+2,mj) =           sum3(i1)
            end do
    
            if (j <= this%jmax) then
              ijml = 3*(j*(j+1)/2+m)
    
              if (j == 0) then
                cc(19,mj) = -v(ijml+1)
              else
                cc(19,mj) = sqrt(j/(2*j+1._dbl)) * v(ijml-1) - sqrt((j+1)/(2*j+1._dbl)) * v(ijml+1)
              end if
            end if
        end do
      end do
      
      do mj = 1, this%jms2
        do i1 = 1, 5
          cc(i1,mj) = cc(i1,mj) / 2
        end do
    
        do i1 = 7, 18
          cc(i1,mj) = cc(i1,mj) / 2
        end do
      end do
    
    deallocate(cab, sum1, sum2, sum3)
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), fftLege(step), sinx(step), symL(step,19), asymL(step,19), &
            & sumLegendreN(19,step,0:this%maxj), sumLegendreS(19,step,0:this%maxj), grid(19,step,0:this%nFourier-1),  &
            & cr(4,this%jms1), fft(4,step,0:this%nFourier-1), sumFourierN(4,step,0:this%jmax+1), symF(step,4), asymF(step,4), &
            & sumFourierS(4,step,0:this%jmax+1) )
      
      cr = czero

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
            
            do concurrent ( i1=1:19 , i2=1:step )
              symL(i2,i1) = symL(i2,i1) + cc(i1,mj)
            end do
            
          do j = 1, (this%maxj-m)/2
            mj = mj+2
            
            pmj2 = pmj1
            pmj1 = pmj
            pmj = this%amjrr(mj-1) * cosx * pmj1 - this%bmjrr(mj-1) * pmj2
            
            do concurrent ( i1=1:19 , i2=1:step )
              asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj-1) * pmj(i2)
            end do
            
            pmj2 = pmj1
            pmj1 = pmj
            pmj = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
            
            do concurrent ( i1=1:19 , i2=1:step )
              symL(i2,i1) = symL(i2,i1) + cc(i1,mj) * pmj(i2)
            end do
            
            if ( maxval(abs(pmj)) < this%tolm ) exit
          end do
          
          if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
            mj = mj+1
            
            pmj2 = pmj1
            pmj1 = pmj
            pmj = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
            
            do concurrent ( i1=1:19 , i2=1:step )
              asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj) * pmj(i2)
            end do
          end if
          
          do concurrent ( i2=1:step, i1=1:19 )
            sumLegendreN(i1,i2,0)%re = symL(i2,i1)%re + asymL(i2,i1)%re
            sumLegendreS(i1,i2,0)%re = symL(i2,i1)%re - asymL(i2,i1)%re
            
            sumLegendreN(i1,i2,0)%im = 0._dbl ; sumLegendreS(i1,i2,0)%im = 0._dbl
          end do
          
        do m = 1, this%maxj
          fac = -sqrt( ( 2*m+1 ) / ( 2._dbl*m ) )
          pmm = fac * sinx * pmm ; if (maxval(abs(pmm)) < this%tolm) exit
          
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = pmm
          
          symL  = czero
          asymL = czero

          j = m
            mj = m*(this%maxj+1)-m*(m+1)/2+j+1
            
            do concurrent ( i1=1:19 , i2=1:step )
              symL(i2,i1) = symL(i2,i1) + cc(i1,mj) * pmj(i2)
            end do
          
          do j = 1, (this%maxj-m)/2
            mj = mj+2

            pmj2 = pmj1
            pmj1 = pmj
            pmj = this%amjrr(mj-1) * cosx * pmj1 - this%bmjrr(mj-1) * pmj2

            do concurrent ( i1=1:19 , i2=1:step )
              asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj-1) * pmj(i2)
            end do

            pmj2 = pmj1
            pmj1 = pmj
            pmj  = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2

            do concurrent ( i1=1:19 , i2=1:step )
              symL(i2,i1) = symL(i2,i1) + cc(i1,mj) * pmj(i2)
            end do
            
            if ( maxval(abs(pmj)) < this%tolm ) exit
          end do
          
          if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
            mj = mj+1

            pmj2 = pmj1
            pmj1 = pmj
            pmj  = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2

            do concurrent ( i1=1:19 , i2=1:step )
              asymL(i2,i1) = asymL(i2,i1) + cc(i1,mj) * pmj(i2)
            end do
          end if
          
          do concurrent ( i2=1:step, i1=1:19 )
            sumLegendreN(i1,i2,m) = ( symL(i2,i1) + asymL(i2,i1) )
            sumLegendreS(i1,i2,m) = ( symL(i2,i1) - asymL(i2,i1) ) 
          end do
        end do
        
        call this%fourtrans%exec_c2r_sub(19*step, this%maxj+1, sumLegendreN, grid)
        
        do concurrent (i1=0:this%nFourier-1, i2=1:step)
          fft(1,i2,i1) = grid( 1,i2,i1) * grid( 4,i2,i1) + &    
                       & grid( 2,i2,i1) * grid( 5,i2,i1) + &    
                       & grid( 3,i2,i1) * grid( 6,i2,i1)
          fft(2,i2,i1) = grid( 4,i2,i1) * grid( 7,i2,i1) + &
                       & grid( 5,i2,i1) * grid( 8,i2,i1) + &
                       & grid( 6,i2,i1) * grid( 9,i2,i1) + &
                       & grid(16,i2,i1) * grid(19,i2,i1)
          fft(3,i2,i1) = grid( 4,i2,i1) * grid(10,i2,i1) + &
                       & grid( 5,i2,i1) * grid(11,i2,i1) + &
                       & grid( 6,i2,i1) * grid(12,i2,i1) + &
                       & grid(17,i2,i1) * grid(19,i2,i1)
          fft(4,i2,i1) = grid( 4,i2,i1) * grid(13,i2,i1) + &
                       & grid( 5,i2,i1) * grid(14,i2,i1) + &
                       & grid( 6,i2,i1) * grid(15,i2,i1) + &
                       & grid(18,i2,i1) * grid(19,i2,i1)
        end do

        call this%fourtrans%exec_r2c_sub(4*step, this%maxj, fft, sumFourierN) ; sumFourierN = sumFourierN * this%nFourier

        call this%fourtrans%exec_c2r_sub(19*step, this%maxj+1, sumLegendreS, grid)
        
        do concurrent (i1=0:this%nFourier-1, i2=1:step)
          fft(1,i2,i1) = grid( 1,i2,i1) * grid( 4,i2,i1) + &    
                       & grid( 2,i2,i1) * grid( 5,i2,i1) + &    
                       & grid( 3,i2,i1) * grid( 6,i2,i1)
          fft(2,i2,i1) = grid( 4,i2,i1) * grid( 7,i2,i1) + &
                       & grid( 5,i2,i1) * grid( 8,i2,i1) + &
                       & grid( 6,i2,i1) * grid( 9,i2,i1) + &
                       & grid(16,i2,i1) * grid(19,i2,i1)
          fft(3,i2,i1) = grid( 4,i2,i1) * grid(10,i2,i1) + &
                       & grid( 5,i2,i1) * grid(11,i2,i1) + &
                       & grid( 6,i2,i1) * grid(12,i2,i1) + &
                       & grid(17,i2,i1) * grid(19,i2,i1)
          fft(4,i2,i1) = grid( 4,i2,i1) * grid(13,i2,i1) + &
                       & grid( 5,i2,i1) * grid(14,i2,i1) + &
                       & grid( 6,i2,i1) * grid(15,i2,i1) + &
                       & grid(18,i2,i1) * grid(19,i2,i1)
        end do
        
        call this%fourtrans%exec_r2c_sub(4*step, this%maxj, fft, sumFourierS) ; sumFourierS = sumFourierS * this%nFourier
        
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
          fac = -sqrt( ( 2*m+1 ) / ( 2._dbl*m ) )
          pmm = fac * sinx * pmm
          
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
    
    fac = 1 / ( 4 * this%nLegendre**2 * this%nFourier * sqrt(pi) )
    
    do mj = 1, this%jms1 
      cr(1,mj) = cr(1,mj) * fac
      cr(2,mj) = cr(2,mj) * fac / 2
      cr(3,mj) = cr(3,mj) * fac / 2 * cunit
      cr(4,mj) = cr(4,mj) * fac

      cr12     = cr(2,mj) - cr(3,mj)
      cr(3,mj) = cr(2,mj) + cr(3,mj)
      cr(2,mj) = cr12
    end do
    
    j = 0
      m = 0
        ijm = 1
        lm  = m*(this%maxj)-m*(m+1)/2+j+1
        lm2 = lm+this%maxj-m

        cjm(1,ijm) =        cr(1,lm  )                                ; cjm(1,ijm)%im = 0._dbl
        cjm(4,ijm) =        cr(2,lm2 )   * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                   &        cr(4,lm+1)   * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                   & conjg( cr(2,lm2 ) ) * cleb1_fn(j+1,m-1,1,+1,j,m) ; cjm(4,ijm)%im = 0._dbl

    do j = 1, this%jmax
      m = 0
        ijm = ijm+1
        lm  = m*(this%maxj)-m*(m+1)/2+j+1
        lm2 = lm+this%maxj-m-1

        cjm(1,ijm) =        cr(1,lm   )                                ; cjm(1,ijm)%im = 0._dbl
        cjm(2,ijm) =        cr(2,lm2-1)   * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                   &        cr(4,lm -1)   * cleb1_fn(j-1,m+0,1, 0,j,m) + &
                   & conjg( cr(2,lm2-1) ) * cleb1_fn(j-1,m-1,1,+1,j,m) ; cjm(2,ijm)%im = 0._dbl
        cjm(3,ijm) =        cr(2,lm2  )   * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                   &        cr(4,lm   )   * cleb1_fn(j  ,m+0,1, 0,j,m) + &
                   & conjg( cr(2,lm2  ) ) * cleb1_fn(j  ,m-1,1,+1,j,m) ; cjm(3,ijm)%re = 0._dbl
        cjm(4,ijm) =        cr(2,lm2+1)   * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                   &        cr(4,lm +1)   * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                   & conjg( cr(2,lm2+1) ) * cleb1_fn(j+1,m-1,1,+1,j,m) ; cjm(4,ijm)%im = 0._dbl

      do m = 1, j
        ijm = ijm+1
        lm  = m*(this%maxj)-m*(m+1)/2+j+1
        lm1 = lm - this%maxj + m
        lm2 = lm + this%maxj - m - 1

        cjm(1,ijm) = cr(1,lm   )
        cjm(2,ijm) = cr(2,lm2-1) * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                   & cr(4,lm -1) * cleb1_fn(j-1,m+0,1, 0,j,m) - &
                   & cr(3,lm1-1) * cleb1_fn(j-1,m-1,1,+1,j,m)
        cjm(3,ijm) = cr(2,lm2  ) * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                   & cr(4,lm   ) * cleb1_fn(j  ,m+0,1, 0,j,m) - &
                   & cr(3,lm1  ) * cleb1_fn(j  ,m-1,1,+1,j,m)
        cjm(4,ijm) = cr(2,lm2+1) * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                   & cr(4,lm +1) * cleb1_fn(j+1,m+0,1, 0,j,m) - &
                   & cr(3,lm1+1) * cleb1_fn(j+1,m-1,1,+1,j,m)
      end do
    end do

    deallocate(cr)
    
  end subroutine vcsv_vcvv_vcvgv_sub

  subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub(this)
    class(T_lateralGrid), intent(inout) :: this

    call this%fourtrans%deallocate_sub()

  end subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub
  
end submodule vcsv_vcvv_vcvgv
