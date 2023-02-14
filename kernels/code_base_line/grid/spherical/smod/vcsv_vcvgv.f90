submodule (SphericalHarmonics) vcsv_vcvgv
  implicit none

  contains

  subroutine init_vcsv_vcvgv_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer,              allocatable   :: iemb(:), oemb(:)
    real(kind=dbl),       allocatable   :: testField_re(:,:), testField3_re(:,:,:)
    complex(kind=dbl),    allocatable   :: testField(:,:), testField3(:,:,:)

    allocate(testField3(16,step,this%nFourier))
      this%fftw_16_forw = fftw_plan_many_dft( 1, (/this%nFourier/), 16*step, testField3, iemb, 16*step, 1,                &
                                            &                                testField3, oemb, 16*step, 1, +1, fftw_flags )
    deallocate(testField3)

    allocate(testField3_re(step,3,this%nFourier), testField3(step,3,this%nFourier/2+1))
      this%fftw_03_back = fftw_plan_many_dft_r2c( 1, (/this%nFourier/), 3*step, testField3_re, iemb, 3*step, 1,            &
                                                &                               testField3   , oemb, 3*step, 1, fftw_flags )
    deallocate(testField3_re, testField3)

    write(*,*) 'vcsv_vcvgv initialized'

  end subroutine init_vcsv_vcvgv_sub

  function vcsv_vcvgv_fn(this, ri, dv_r, v) result(cjml)
    class(T_lateralGrid), intent(in) :: this
    real(kind=dbl),       intent(in) :: ri
    complex(kind=dbl),    intent(in) :: v(:), dv_r(:)
    complex(kind=dbl)                :: cjml(this%jmv)
    integer                          :: i, k, j, m, l, jm_int, lm_int, jml_int, i1, i2, mj
    real(kind=dbl),     allocatable  :: gc(:), pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), fftLege(:), fft(:,:,:)
    complex(kind=dbl),  allocatable  :: sum1(:), sum2(:), sum3(:), cab(:,:), cc(:,:), cr(:,:)
    complex(kind=dbl),  allocatable  :: sumLegendreN(:,:,:), sumLegendreS(:,:,:), fftNC(:,:,:), fftSC(:,:,:)

    allocate(cab(5,this%jmv1), sum1(2), sum2(2), sum3(2), gc(2)) ; cab = cmplx(0._dbl, 0._dbl, kind=dbl)
      
      do jml_int = 1, this%jmv
        cab(1,jml_int) = v(jml_int)
        cab(5,jml_int) = dv_r(jml_int)
      end do

      do j = 1, this%jmax+1
        gc(1) = sqrt(j*(j+1._dbl)*(j+1._dbl)/(2*j+1._dbl))
        gc(2) = sqrt(j*(j       )*(j+1._dbl)/(2*j+1._dbl))

        do m = 0, j
          sum1 = cmplx(0._dbl, 0._dbl, kind=dbl)
          sum2 = cmplx(0._dbl, 0._dbl, kind=dbl)
          sum3 = cmplx(0._dbl, 0._dbl, kind=dbl)

          if (m == 0) then
            do l = abs(j-1), min(this%jmax, j+1)
              sum1 = sum1 + conjg( v(3*(l*(l+1)/2+m+1)+j-l) ) * cleb1_fn(j,m,1,-1,l,m-1) * gc * (-1)**(j+l)
              sum3 = sum3 +        v(3*(l*(l+1)/2+m+0)+j-l)   * cleb1_fn(j,m,1, 0,l,m+0) * gc
              sum2 = sum2 +        v(3*(l*(l+1)/2+m+1)+j-l)   * cleb1_fn(j,m,1,+1,l,m+1) * gc
            end do
          else
            do l = max(abs(m-1), j-1), min(this%jmax, j+1)
                           sum1 = sum1 + v(3*(l*(l+1)/2+m-1)+j-l) * cleb1_fn(j,m,1,-1,l,m-1) * gc
              if (l > m-1) sum3 = sum3 + v(3*(l*(l+1)/2+m+0)+j-l) * cleb1_fn(j,m,1, 0,l,m+0) * gc
              if (l > m+0) sum2 = sum2 + v(3*(l*(l+1)/2+m+1)+j-l) * cleb1_fn(j,m,1,+1,l,m+1) * gc
            end do
          end if

          jml_int = 3*(j*(j+1)/2+m)+abs(j-1)-j
            cab(2,jml_int) =         sum1(1) - sum2(1)
            cab(3,jml_int) = cunit*( sum1(1) + sum2(1) )
            cab(4,jml_int) =         sum3(1)

          jml_int = 3*(j*(j+1)/2+m)+(j+1)-j
            cab(2,jml_int) =         sum1(2) - sum2(2)
            cab(3,jml_int) = cunit*( sum1(2) + sum2(2) )
            cab(4,jml_int) =         sum3(2)
        end do
      end do

      cab(2:4,:) = cab(2:4,:) / ri

    deallocate(sum1, sum2, sum3, gc)
    allocate(cc(16, this%jms2), sum1(5), sum2(5), sum3(5)); cc = cmplx(0._dbl, 0._dbl, kind=dbl)
    
      mj = 0
        do m = 0, this%maxj
          do j = m, this%maxj
            sum1 = cmplx(0._dbl, 0._dbl, kind=dbl)
            sum2 = cmplx(0._dbl, 0._dbl, kind=dbl)
            sum3 = cmplx(0._dbl, 0._dbl, kind=dbl)

            do l = abs(j-1), min(this%jmax+1, j+1)
              if ( m == 0 ) then
                sum1 = sum1 + conjg( cab(:,3*(l*(l+1)/2+m+1)+j-l) ) * cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(j+l)
                sum3 = sum3 +        cab(:,3*(l*(l+1)/2+m  )+j-l)   * cleb1_fn(j,m,1, 0,l,m  )
                sum2 = sum2 +        cab(:,3*(l*(l+1)/2+m+1)+j-l)   * cleb1_fn(j,m,1,+1,l,m+1)
              else
                             sum1 = sum1 + cab(:,3*(l*(l+1)/2+m-1)+j-l) * cleb1_fn(j,m,1,-1,l,m-1)
                if (l > m-1) sum3 = sum3 + cab(:,3*(l*(l+1)/2+m  )+j-l) * cleb1_fn(j,m,1, 0,l,m  )
                if (l > m+0) sum2 = sum2 + cab(:,3*(l*(l+1)/2+m+1)+j-l) * cleb1_fn(j,m,1,+1,l,m+1)
              end if
            end do

            mj = mj+1
              do i1 = 1, 5
                cc(3*(i1-1)+1, mj) =           sum1(i1) - sum2(i1)
                cc(3*(i1-1)+2, mj) = cunit * ( sum1(i1) + sum2(i1) )
                cc(3*(i1-1)+3, mj) =           sum3(i1)
              end do
              
              if (j <= this%jmax) then
                jml_int = 3*(j*(j+1)/2+m)
  
                if (j == 0) then
                  cc(16,mj) = -v(jml_int+1)
                else
                  cc(16,mj) = sqrt(j/(2*j+1._dbl)) * v(jml_int-1) - sqrt((j+1)/(2*j+1._dbl)) * v(jml_int+1)
                end if
              end if
          end do
        end do
      
      cc(:,1:this%maxj+1) = cc(:,1:this%maxj+1) / 2 ; cc(3,:) = cc(3,:) * 2 ; cc(16,:) = cc(16,:) * 2

    deallocate(cab, sum1, sum2, sum3)
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), sinx(step), fftLege(step),             &
            & sumLegendreN(16,step,0:this%nFourier-1), sumLegendreS(16,step,0:this%nFourier-1),                &
            & fft(step,3,0:this%nFourier-1), fftNC(step,3,0:this%nFourier/2), fftSC(step,3,0:this%nFourier/2), &
            & cr(3,this%jms2) ) ; cr = cmplx(0._dbl, 0._dbl, kind=dbl)

      do i = 1, this%nLegendre, step
        cosx    = this%roots(i:i+step-1)
        sinx    = sqrt(1-cosx**2)
        fftLege = this%fftLege(i:i+step-1)

        sumLegendreN = cmplx(0._dbl, 0._dbl, kind=dbl)
        sumLegendreS = cmplx(0._dbl, 0._dbl, kind=dbl)
        
        pmm = 1._dbl; mj = 0
          do m = 0, this%maxj
            pmj2 = 0._dbl
            pmj1 = 0._dbl
            pmj  = 1._dbl

            do j = m, this%maxj
              mj = mj+1

              do concurrent (i2=1:step, i1=1:16)
                sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj) * pmj(i2)
                sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj) * pmj(i2) * (-1)**(j+m)
              end do

              pmj2 = pmj1
              pmj1 = pmj
              pmj  = this%amjrr(mj+1) * cosx * pmj1 - this%bmjrr(mj+1) * pmj2
            end do

            do concurrent (i2=1:step, i1=1:16)
              sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) * pmm(i2)
              sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) * pmm(i2)
            end do

            pmm = -sinx * this%cmmrr(m) * pmm; if (maxval(abs(pmm)) < 1.0d-55) exit
          end do

        call fftw_execute_dft(this%fftw_16_forw, sumLegendreN, sumLegendreN)
          do concurrent (i1 = 0:this%nFourier-1, i2 = 1:step)
            fft(i2,1,i1) = sumLegendreN( 1,i2,i1)%re * sumLegendreN( 4,i2,i1)%re + &
                         & sumLegendreN( 2,i2,i1)%re * sumLegendreN( 5,i2,i1)%re + &
                         & sumLegendreN( 3,i2,i1)%re * sumLegendreN( 6,i2,i1)%re + &
                         & sumLegendreN(13,i2,i1)%re * sumLegendreN(16,i2,i1)%re
            fft(i2,2,i1) = sumLegendreN( 1,i2,i1)%re * sumLegendreN( 7,i2,i1)%re + &
                         & sumLegendreN( 2,i2,i1)%re * sumLegendreN( 8,i2,i1)%re + &
                         & sumLegendreN( 3,i2,i1)%re * sumLegendreN( 9,i2,i1)%re + &
                         & sumLegendreN(14,i2,i1)%re * sumLegendreN(16,i2,i1)%re
            fft(i2,3,i1) = sumLegendreN( 1,i2,i1)%re * sumLegendreN(10,i2,i1)%re + &
                         & sumLegendreN( 2,i2,i1)%re * sumLegendreN(11,i2,i1)%re + &
                         & sumLegendreN( 3,i2,i1)%re * sumLegendreN(12,i2,i1)%re + &
                         & sumLegendreN(15,i2,i1)%re * sumLegendreN(16,i2,i1)%re
          end do
        call fftw_execute_dft_r2c(this%fftw_03_back, fft, fftNC)

        call fftw_execute_dft(this%fftw_16_forw, sumLegendreS, sumLegendreS)
          do concurrent (i1 = 0:this%nFourier-1, i2 = 1:step)
            fft(i2,1,i1) = sumLegendreS( 1,i2,i1)%re * sumLegendreS( 4,i2,i1)%re + &
                         & sumLegendreS( 2,i2,i1)%re * sumLegendreS( 5,i2,i1)%re + &
                         & sumLegendreS( 3,i2,i1)%re * sumLegendreS( 6,i2,i1)%re + &
                         & sumLegendreS(13,i2,i1)%re * sumLegendreS(16,i2,i1)%re
            fft(i2,2,i1) = sumLegendreS( 1,i2,i1)%re * sumLegendreS( 7,i2,i1)%re + &
                         & sumLegendreS( 2,i2,i1)%re * sumLegendreS( 8,i2,i1)%re + &
                         & sumLegendreS( 3,i2,i1)%re * sumLegendreS( 9,i2,i1)%re + &
                         & sumLegendreS(14,i2,i1)%re * sumLegendreS(16,i2,i1)%re
            fft(i2,3,i1) = sumLegendreS( 1,i2,i1)%re * sumLegendreS(10,i2,i1)%re + &
                         & sumLegendreS( 2,i2,i1)%re * sumLegendreS(11,i2,i1)%re + &
                         & sumLegendreS( 3,i2,i1)%re * sumLegendreS(12,i2,i1)%re + &
                         & sumLegendreS(15,i2,i1)%re * sumLegendreS(16,i2,i1)%re
          end do
        call fftw_execute_dft_r2c(this%fftw_03_back, fft, fftSC)

        pmm = 1._dbl; mj = 1
          do m = 0, this%jmax+1
            pmj2 = 0._dbl
            pmj1 = 0._dbl
            pmj  = 1._dbl

            do concurrent ( i1 = 1:3, i2 = 1:step )
              fftNC(i2,i1,m) = fftLege(i2) * fftNC(i2,i1,m) * pmm(i2)
              fftSC(i2,i1,m) = fftLege(i2) * fftSC(i2,i1,m) * pmm(i2)
            end do

            do j = m, this%jmax+1
              jm_int = j*(j+1)/2+m+1
                do i1 = 1, 3
                  cr(i1,jm_int) = cr(i1,jm_int) + sum( pmj(:) * (fftNC(:,i1,m) + (-1)**(j+m) * fftSC(:,i1,m)) )
                end do

              mj = mj+1
                pmj2 = pmj1
                pmj1 = pmj
                pmj  = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
            end do
          
            pmm = -sinx * this%cmmrr(m) * pmm ; if (maxval(abs(pmm)) < 1.0d-55) exit
          end do
        end do

        cr = cr / 4 / this%nLegendre**2 / this%nFourier / sqrt(pi)

    deallocate( cc, fft, fftNC, fftSC, sumLegendreN, sumLegendreS, pmm, pmj, pmj1, pmj2, cosx, sinx, fftLege )
    
      do j = 0, this%jmax
        do m = 0, j
          do l = abs(j-1), j+1
            lm_int = l*(l+1)/2+m+1; jml_int = 3*(j*(j+1)/2+m)+l-j
            if ( m == 0 ) then
              cjml(jml_int) = (       cr(1,lm_int+1) - cunit * cr(2,lm_int+1)  ) * cleb1_fn(l,m+1,1,-1,j,m) / 2 + &
                            & (       cr(3,lm_int  )                           ) * cleb1_fn(l,m+0,1, 0,j,m)     + &
                            & ( conjg(cr(1,lm_int+1) - cunit * cr(2,lm_int+1)) ) * cleb1_fn(l,m-1,1,+1,j,m) / 2
            else
              cjml(jml_int) = ( cr(1,lm_int+1) - cunit * cr(2,lm_int+1) ) * cleb1_fn(l,m+1,1,-1,j,m) / 2 + &
                            & ( cr(3,lm_int  )                          ) * cleb1_fn(l,m+0,1, 0,j,m)     - &
                            & ( cr(1,lm_int-1) + cunit * cr(2,lm_int-1) ) * cleb1_fn(l,m-1,1,+1,j,m) / 2
            end if
          end do
        end do
      end do

    deallocate(cr)

  end function vcsv_vcvgv_fn

end submodule vcsv_vcvgv