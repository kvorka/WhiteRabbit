submodule (SphericalHarmonics) vcsv_vcvv_vcvgv
  implicit none

  contains

  subroutine init_vcsv_vcvv_vcvgv_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer,              pointer       :: ip(:) => null()
    complex(kind=dbl),    allocatable   :: field(:,:,:)

    allocate(field(19,step,this%nFourier))
      this%fftw_19_forw = fftw_plan_many_dft( 1, [this%nFourier], 19*step, field, ip, 19*step, 1,                &
                                            &                              field, ip, 19*step, 1, +1, fftw_flags )
    deallocate(field)

    allocate(field(8,step,this%nFourier/2))
      this%fftw_08_back = fftw_plan_many_dft( 1, [this%nFourier]/2, 8*step, field, ip, 8*step, 1,                &
                                            &                               field, ip, 8*step, 1, -1, fftw_flags )
    deallocate(field)

    write(*,*) 'vcsv_vcvv_vcvgv initialized'

  end subroutine init_vcsv_vcvv_vcvgv_sub

  subroutine vcsv_vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm, cjml)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(:), q(:), v(:)
    complex(kind=dbl),    intent(out) :: cjm(:), cjml(:)
    integer                           :: i, j, m, l, jm_int, lm, jml_int, i1, i2, mj, s, iL, lm1, lm2
    real(kind=dbl)                    :: cmm
    complex(kind=dbl)                 :: mult, jexp
    real(kind=dbl),       allocatable :: gc(:), pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), fftLege(:)
    complex(kind=dbl),    allocatable :: sum1(:), sum2(:), sum3(:), cab(:,:), cc(:,:), cr(:,:), symL(:,:), asymL(:,:)
    complex(kind=dbl),    allocatable :: sumLegendreN(:,:,:), sumLegendreS(:,:,:), fftC(:,:), fft(:,:,:)

    allocate(cab(6,this%jmv1), sum1(2), sum2(2), sum3(2), gc(2)) ; cab = cmplx(0._dbl, 0._dbl, kind=dbl)

      do jml_int = 1, this%jmv
        cab(1,jml_int) = q(jml_int)
        cab(2,jml_int) = v(jml_int)
        cab(6,jml_int) = dv_r(jml_int)
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
            cab(3,jml_int) =        sum1(1) - sum2(1)
            cab(4,jml_int) = cunit*(sum1(1) + sum2(1))
            cab(5,jml_int) =        sum3(1)

          jml_int = 3*(j*(j+1)/2+m)+(j+1)-j
            cab(3,jml_int) =        sum1(2) - sum2(2)
            cab(4,jml_int) = cunit*(sum1(2) + sum2(2))
            cab(5,jml_int) =        sum3(2)
        end do
      end do

      cab(3:5,:) = cab(3:5,:) / ri

    deallocate(sum1, sum2, sum3, gc)
    allocate(cc(19, this%jms2), sum1(6), sum2(6), sum3(6)) ; cc = cmplx(0._dbl, 0._dbl, kind=dbl)
    
    mj = 0
      do m = 0, this%maxj
        do j = m, this%maxj
          sum1 = cmplx(0._dbl, 0._dbl, kind=dbl)
          sum2 = cmplx(0._dbl, 0._dbl, kind=dbl)
          sum3 = cmplx(0._dbl, 0._dbl, kind=dbl)

          if (m == 0) then
            do l = abs(j-1), min(this%jmax+1, j+1)
              sum1 = sum1 + conjg( cab(:,3*(l*(l+1)/2+m+1)+j-l) ) * cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(j+l)
              sum3 = sum3 +        cab(:,3*(l*(l+1)/2+m  )+j-l)   * cleb1_fn(j,m,1, 0,l,m  )
              sum2 = sum2 +        cab(:,3*(l*(l+1)/2+m+1)+j-l)   * cleb1_fn(j,m,1,+1,l,m+1)
            end do
          else
            do l = abs(j-1), min(this%jmax+1, j+1)
                           sum1 = sum1 + cab(:,3*(l*(l+1)/2+m-1)+j-l) * cleb1_fn(j,m,1,-1,l,m-1)
              if (l > m-1) sum3 = sum3 + cab(:,3*(l*(l+1)/2+m  )+j-l) * cleb1_fn(j,m,1, 0,l,m  )
              if (l > m+0) sum2 = sum2 + cab(:,3*(l*(l+1)/2+m+1)+j-l) * cleb1_fn(j,m,1,+1,l,m+1)
            end do
          end if

          mj = mj+1
            do i1 = 1, 6
              cc(3*(i1-1)+1, mj) =           sum1(i1) - sum2(i1)
              cc(3*(i1-1)+2, mj) = cunit * ( sum1(i1) + sum2(i1) )
              cc(3*(i1-1)+3, mj) =           sum3(i1)
            end do

            if (j <= this%jmax) then
              jml_int = 3*(j*(j+1)/2+m)

              if (j == 0) then
                cc(19,mj) = -v(jml_int+1)
              else
                cc(19,mj) = sqrt(j/(2*j+1._dbl)) * v(jml_int-1) - sqrt((j+1)/(2*j+1._dbl)) * v(jml_int+1)
              end if
            end if
        end do
      end do
      
      cc(:,1:this%maxj+1) = cc(:,1:this%maxj+1) / 2 ; cc(6,:) = 2 * cc(6,:) ; cc(19,:) = 2 * cc(19,:)

    deallocate(cab, sum1, sum2, sum3)
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), sinx(step), fftLege(step),                             &
            & symL(19,step), asymL(19,step), sumLegendreN(19,step,0:this%nFourier-1), sumLegendreS(19,step,0:this%nFourier-1), &                                 
            & fft(8,step,0:this%nFourier/2-1), fftC(8,step), cr(4,this%jms1)                                                   )
      
      cr = cmplx(0._dbl, 0._dbl, kind=dbl)

      do i = 1, this%nLegendre, step
        cosx    = this%roots(i:i+step-1)
        sinx    = sqrt(1 - cosx**2)
        fftLege = this%fftLege(i:i+step-1)
        
        pmm = 1._dbl
        mj  = 0

        sumLegendreN = cmplx(0._dbl, 0._dbl, kind=dbl)
        sumLegendreS = cmplx(0._dbl, 0._dbl, kind=dbl)

        m = 0
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl
          
          symL  = cmplx(0._dbl, 0._dbl, kind=dbl)
          asymL = cmplx(0._dbl, 0._dbl, kind=dbl)
          
          j = m
            s = +1 ; mj = mj+1 
            
            do concurrent ( i2=1:step, i1=1:19 )
              symL(i1,i2) = symL(i1,i2) + cc(i1,mj)
            end do
            
          do j = m+1, this%maxj
            s = -s ; mj = mj+1

            pmj2 = pmj1
            pmj1 = pmj
            pmj  = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
            
            if (s == 1) then
              do concurrent ( i2=1:step, i1=1:19 )
                symL(i1,i2) = symL(i1,i2) + cc(i1,mj) * pmj(i2)
              end do
            else
              do concurrent ( i2=1:step, i1=1:19 )
                asymL(i1,i2) = asymL(i1,i2) + cc(i1,mj) * pmj(i2)
              end do
            end if
          end do
          
          do concurrent ( i2=1:step, i1=1:19 )
            sumLegendreN(i1,i2,0)%re = symL(i1,i2)%re+asymL(i1,i2)%re ; sumLegendreN(i1,i2,0)%im = 0._dbl
            sumLegendreS(i1,i2,0)%re = symL(i1,i2)%re-asymL(i1,i2)%re ; sumLegendreS(i1,i2,0)%im = 0._dbl
          end do
          
        do m = 1, this%maxj
          cmm = -sqrt( ( 2*m+1 ) / ( 2._dbl*m ) )
          pmm = cmm * sinx * pmm
          if (maxval(abs(pmm)) < 1.0d-55) exit
          
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl
          
          symL  = cmplx(0._dbl, 0._dbl, kind=dbl)
          asymL = cmplx(0._dbl, 0._dbl, kind=dbl)

          j = m
          s = +1 ; mj = mj+1
            
            do concurrent ( i2=1:step, i1=1:19 )
              symL(i1,i2) = symL(i1,i2) + cc(i1,mj)
            end do
          
          do j = m+1, this%maxj
            s = -s ; mj = mj+1

            pmj2 = pmj1
            pmj1 = pmj
            pmj  = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2
            
            if (s == 1) then
              do concurrent ( i2=1:step, i1=1:19 )
                symL(i1,i2) = symL(i1,i2) + cc(i1,mj) * pmj(i2)
              end do
            else
              do concurrent ( i2=1:step, i1=1:19 )
                asymL(i1,i2) = asymL(i1,i2) + cc(i1,mj) * pmj(i2)
              end do
            end if
          end do
          
          do concurrent ( i2=1:step, i1=1:19 )
            sumLegendreN(i1,i2,m) = (symL(i1,i2)+asymL(i1,i2)) * pmm(i2)
            sumLegendreS(i1,i2,m) = (symL(i1,i2)-asymL(i1,i2)) * pmm(i2)
          end do
        end do

        call fftw_execute_dft(this%fftw_19_forw, sumLegendreN, sumLegendreN)
        call fftw_execute_dft(this%fftw_19_forw, sumLegendreS, sumLegendreS)
        
        do concurrent (i1=0:this%nFourier/2-1, i2=1:step)
          fft(1,i2,i1)%re = sumLegendreN( 1,i2,2*i1  )%re * sumLegendreN( 4,i2,2*i1  )%re + &    
                          & sumLegendreN( 2,i2,2*i1  )%re * sumLegendreN( 5,i2,2*i1  )%re + &    
                          & sumLegendreN( 3,i2,2*i1  )%re * sumLegendreN( 6,i2,2*i1  )%re
          fft(3,i2,i1)%re = sumLegendreN( 4,i2,2*i1  )%re * sumLegendreN( 7,i2,2*i1  )%re + &
                          & sumLegendreN( 5,i2,2*i1  )%re * sumLegendreN( 8,i2,2*i1  )%re + &
                          & sumLegendreN( 6,i2,2*i1  )%re * sumLegendreN( 9,i2,2*i1  )%re + &
                          & sumLegendreN(16,i2,2*i1  )%re * sumLegendreN(19,i2,2*i1  )%re
          fft(5,i2,i1)%re = sumLegendreN( 4,i2,2*i1  )%re * sumLegendreN(10,i2,2*i1  )%re + &
                          & sumLegendreN( 5,i2,2*i1  )%re * sumLegendreN(11,i2,2*i1  )%re + &
                          & sumLegendreN( 6,i2,2*i1  )%re * sumLegendreN(12,i2,2*i1  )%re + &
                          & sumLegendreN(17,i2,2*i1  )%re * sumLegendreN(19,i2,2*i1  )%re
          fft(7,i2,i1)%re = sumLegendreN( 4,i2,2*i1  )%re * sumLegendreN(13,i2,2*i1  )%re + &
                          & sumLegendreN( 5,i2,2*i1  )%re * sumLegendreN(14,i2,2*i1  )%re + &
                          & sumLegendreN( 6,i2,2*i1  )%re * sumLegendreN(15,i2,2*i1  )%re + &
                          & sumLegendreN(18,i2,2*i1  )%re * sumLegendreN(19,i2,2*i1  )%re
          
          fft(1,i2,i1)%im = sumLegendreN( 1,i2,2*i1+1)%re * sumLegendreN( 4,i2,2*i1+1)%re + &    
                          & sumLegendreN( 2,i2,2*i1+1)%re * sumLegendreN( 5,i2,2*i1+1)%re + &    
                          & sumLegendreN( 3,i2,2*i1+1)%re * sumLegendreN( 6,i2,2*i1+1)%re
          fft(3,i2,i1)%im = sumLegendreN( 4,i2,2*i1+1)%re * sumLegendreN( 7,i2,2*i1+1)%re + &
                          & sumLegendreN( 5,i2,2*i1+1)%re * sumLegendreN( 8,i2,2*i1+1)%re + &
                          & sumLegendreN( 6,i2,2*i1+1)%re * sumLegendreN( 9,i2,2*i1+1)%re + &
                          & sumLegendreN(16,i2,2*i1+1)%re * sumLegendreN(19,i2,2*i1+1)%re
          fft(5,i2,i1)%im = sumLegendreN( 4,i2,2*i1+1)%re * sumLegendreN(10,i2,2*i1+1)%re + &
                          & sumLegendreN( 5,i2,2*i1+1)%re * sumLegendreN(11,i2,2*i1+1)%re + &
                          & sumLegendreN( 6,i2,2*i1+1)%re * sumLegendreN(12,i2,2*i1+1)%re + &
                          & sumLegendreN(17,i2,2*i1+1)%re * sumLegendreN(19,i2,2*i1+1)%re
          fft(7,i2,i1)%im = sumLegendreN( 4,i2,2*i1+1)%re * sumLegendreN(13,i2,2*i1+1)%re + &
                          & sumLegendreN( 5,i2,2*i1+1)%re * sumLegendreN(14,i2,2*i1+1)%re + &
                          & sumLegendreN( 6,i2,2*i1+1)%re * sumLegendreN(15,i2,2*i1+1)%re + &
                          & sumLegendreN(18,i2,2*i1+1)%re * sumLegendreN(19,i2,2*i1+1)%re
          
          fft(2,i2,i1)%re = sumLegendreS( 1,i2,2*i1  )%re * sumLegendreS( 4,i2,2*i1  )%re + &    
                          & sumLegendreS( 2,i2,2*i1  )%re * sumLegendreS( 5,i2,2*i1  )%re + &    
                          & sumLegendreS( 3,i2,2*i1  )%re * sumLegendreS( 6,i2,2*i1  )%re
          fft(4,i2,i1)%re = sumLegendreS( 4,i2,2*i1  )%re * sumLegendreS( 7,i2,2*i1  )%re + &
                          & sumLegendreS( 5,i2,2*i1  )%re * sumLegendreS( 8,i2,2*i1  )%re + &
                          & sumLegendreS( 6,i2,2*i1  )%re * sumLegendreS( 9,i2,2*i1  )%re + &
                          & sumLegendreS(16,i2,2*i1  )%re * sumLegendreS(19,i2,2*i1  )%re
          fft(6,i2,i1)%re = sumLegendreS( 4,i2,2*i1  )%re * sumLegendreS(10,i2,2*i1  )%re + &
                          & sumLegendreS( 5,i2,2*i1  )%re * sumLegendreS(11,i2,2*i1  )%re + &
                          & sumLegendreS( 6,i2,2*i1  )%re * sumLegendreS(12,i2,2*i1  )%re + &
                          & sumLegendreS(17,i2,2*i1  )%re * sumLegendreS(19,i2,2*i1  )%re
          fft(8,i2,i1)%re = sumLegendreS( 4,i2,2*i1  )%re * sumLegendreS(13,i2,2*i1  )%re + &
                          & sumLegendreS( 5,i2,2*i1  )%re * sumLegendreS(14,i2,2*i1  )%re + &
                          & sumLegendreS( 6,i2,2*i1  )%re * sumLegendreS(15,i2,2*i1  )%re + &
                          & sumLegendreS(18,i2,2*i1  )%re * sumLegendreS(19,i2,2*i1  )%re
          
          fft(2,i2,i1)%im = sumLegendreS( 1,i2,2*i1+1)%re * sumLegendreS( 4,i2,2*i1+1)%re + &    
                          & sumLegendreS( 2,i2,2*i1+1)%re * sumLegendreS( 5,i2,2*i1+1)%re + &    
                          & sumLegendreS( 3,i2,2*i1+1)%re * sumLegendreS( 6,i2,2*i1+1)%re
          fft(4,i2,i1)%im = sumLegendreS( 4,i2,2*i1+1)%re * sumLegendreS( 7,i2,2*i1+1)%re + &
                          & sumLegendreS( 5,i2,2*i1+1)%re * sumLegendreS( 8,i2,2*i1+1)%re + &
                          & sumLegendreS( 6,i2,2*i1+1)%re * sumLegendreS( 9,i2,2*i1+1)%re + &
                          & sumLegendreS(16,i2,2*i1+1)%re * sumLegendreS(19,i2,2*i1+1)%re
          fft(6,i2,i1)%im = sumLegendreS( 4,i2,2*i1+1)%re * sumLegendreS(10,i2,2*i1+1)%re + &
                          & sumLegendreS( 5,i2,2*i1+1)%re * sumLegendreS(11,i2,2*i1+1)%re + &
                          & sumLegendreS( 6,i2,2*i1+1)%re * sumLegendreS(12,i2,2*i1+1)%re + &
                          & sumLegendreS(17,i2,2*i1+1)%re * sumLegendreS(19,i2,2*i1+1)%re
          fft(8,i2,i1)%im = sumLegendreS( 4,i2,2*i1+1)%re * sumLegendreS(13,i2,2*i1+1)%re + &
                          & sumLegendreS( 5,i2,2*i1+1)%re * sumLegendreS(14,i2,2*i1+1)%re + &
                          & sumLegendreS( 6,i2,2*i1+1)%re * sumLegendreS(15,i2,2*i1+1)%re + &
                          & sumLegendreS(18,i2,2*i1+1)%re * sumLegendreS(19,i2,2*i1+1)%re
        end do
        call fftw_execute_dft(this%fftw_08_back, fft, fft)
        
        pmm  = 1._dbl
        mj   = 0
        mult = exp(-2 * pi * cunit / this%nFourier)
        
        m = 0
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl
          
          jexp = cunit
            fftC = ( (1-cunit) * fft(:,:,0) + (1+cunit) * conjg(fft(:,:,0)) ) / 2 ; fftC%im = 0._dbl
          
          do concurrent (i2=1:step, i1=1:8)
            fftC(i1,i2) = fftLege(i2) * fftC(i1,i2)
          end do
          
          j = m
            s = +1 ; mj = mj+1

            do i2 = 1, step
              cr(1,mj) = cr(1,mj) + ( fftC(1,i2) + fftC(2,i2) )
              cr(2,mj) = cr(2,mj) + ( fftC(3,i2) + fftC(4,i2) )
              cr(3,mj) = cr(3,mj) + ( fftC(5,i2) + fftC(6,i2) )
              cr(4,mj) = cr(4,mj) + ( fftC(7,i2) + fftC(8,i2) )
            end do
          
          do j = m+1, this%jmax+1
            s = -s ; mj = mj+1
            
            pmj2 = pmj1
            pmj1 = pmj
            pmj  = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
            
            do i2 = 1, step
              cr(1,mj) = cr(1,mj) + pmj(i2) * ( fftC(1,i2) + s * fftC(2,i2) )
              cr(2,mj) = cr(2,mj) + pmj(i2) * ( fftC(3,i2) + s * fftC(4,i2) )
              cr(3,mj) = cr(3,mj) + pmj(i2) * ( fftC(5,i2) + s * fftC(6,i2) )
              cr(4,mj) = cr(4,mj) + pmj(i2) * ( fftC(7,i2) + s * fftC(8,i2) )
            end do
          end do
        
        do m = 1, this%jmax+1
          cmm = -sqrt( ( 2*m+1 ) / ( 2._dbl*m ) )
          pmm = cmm * sinx * pmm
          if (maxval(abs(pmm)) < 1.0d-55) exit
          
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl
          
          jexp = mult * jexp
            fftC = ( (1-jexp) * fft(:,:,m) + (1+jexp) * conjg(fft(:,:,this%nFourier/2-m)) ) / 2
          
          do concurrent (i2=1:step, i1=1:8)
            fftC(i1,i2) = fftLege(i2) * fftC(i1,i2) * pmm(i2)
          end do
          
          j = m
            s = +1 ; mj = mj+1
            
            do i2 = 1, step
              cr(1,mj) = cr(1,mj) + ( fftC(1,i2) + fftC(2,i2) )
              cr(2,mj) = cr(2,mj) + ( fftC(3,i2) + fftC(4,i2) )
              cr(3,mj) = cr(3,mj) + ( fftC(5,i2) + fftC(6,i2) )
              cr(4,mj) = cr(4,mj) + ( fftC(7,i2) + fftC(8,i2) )
            end do
          
          do j = m+1, this%jmax+1
            s = -s ; mj = mj+1
            
            pmj2 = pmj1
            pmj1 = pmj
            pmj  = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
            
            do i2 = 1, step
              cr(1,mj) = cr(1,mj) + pmj(i2) * ( fftC(1,i2) + s * fftC(2,i2) )
              cr(2,mj) = cr(2,mj) + pmj(i2) * ( fftC(3,i2) + s * fftC(4,i2) )
              cr(3,mj) = cr(3,mj) + pmj(i2) * ( fftC(5,i2) + s * fftC(6,i2) )
              cr(4,mj) = cr(4,mj) + pmj(i2) * ( fftC(7,i2) + s * fftC(8,i2) )
            end do
          end do
        end do
      end do

      cr = cr / 4 / this%nLegendre**2 / this%nFourier / sqrt(pi)

    deallocate( cc, fft, fftC, sumLegendreN, sumLegendreS, pmm, pmj, pmj1, pmj2, cosx, sinx, fftLege, symL, asymL )
    
      do j = 0, this%jmax
        do m = 0, j
          cjm(j*(j+1)/2+m+1) = cr(1,m*(this%maxj)-m*(m+1)/2+j+1)

          do l = abs(j-1), j+1
            jml_int = 3*(j*(j+1)/2+m)+l-j
            lm  = m*(this%maxj)-m*(m+1)/2+l+1
            lm1 = (m-1)*(this%maxj)-(m-1)*(m  )/2+l+1
            lm2 = (m+1)*(this%maxj)-(m+1)*(m+2)/2+l+1

            if ( m == 0 ) then
              cjml(jml_int) = (        cr(2,lm2) - cunit*cr(3,lm2)   ) * cleb1_fn(l,m+1,1,-1,j,m) / 2 + &
                            & (        cr(4,lm )                     ) * cleb1_fn(l,m+0,1, 0,j,m)     + &
                            & ( conjg( cr(2,lm2) - cunit*cr(3,lm2) ) ) * cleb1_fn(l,m-1,1,+1,j,m) / 2
            else
              cjml(jml_int) = ( cr(2,lm2) - cunit*cr(3,lm2) ) * cleb1_fn(l,m+1,1,-1,j,m) / 2 + &
                            & ( cr(4,lm )                   ) * cleb1_fn(l,m+0,1, 0,j,m)     - &
                            & ( cr(2,lm1) + cunit*cr(3,lm1) ) * cleb1_fn(l,m-1,1,+1,j,m) / 2
            end if
          end do
        end do
      end do

    deallocate(cr)

    do j = 0, this%jmax
      jm_int = j*(j+1)/2+1
        cjm(jm_int)%im = 0._dbl
      
      do l = abs(j-1), j+1
        if (j == l) then
          cjml(3*(jm_int-1))%re = 0._dbl
        else
          cjml(3*(jm_int-1)+l-j)%im = 0._dbl
        end if
      end do
    end do

  end subroutine vcsv_vcvv_vcvgv_sub

  subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub(this)
    class(T_lateralGrid), intent(inout) :: this

    call fftw_destroy_plan( this%fftw_19_forw )
    call fftw_destroy_plan( this%fftw_08_back )

  end subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub

end submodule vcsv_vcvv_vcvgv