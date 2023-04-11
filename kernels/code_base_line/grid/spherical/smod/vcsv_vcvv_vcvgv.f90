submodule (SphericalHarmonics) vcsv_vcvv_vcvgv
  implicit none

  contains

  subroutine init_vcsv_vcvv_vcvgv_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer,              pointer       :: ip(:) => null()
    complex(kind=dbl),    allocatable   :: field(:,:,:)

    allocate(field(38,step,this%nFourier))
      this%fftw_38_forw = fftw_plan_many_dft( 1, [this%nFourier], 38*step, field, ip, 38*step, 1,                &
                                            &                              field, ip, 38*step, 1, +1, fftw_flags )
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
    complex(kind=dbl)                 :: mult, jexp
    real(kind=dbl),       allocatable :: gc(:), pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), fftLege(:)
    complex(kind=dbl),    allocatable :: sum1(:), sum2(:), sum3(:), cab(:,:), cc(:,:), cr(:,:), symL(:,:), asymL(:,:)
    complex(kind=dbl),    allocatable :: sumLegendre(:,:,:), fftC(:,:), fft(:,:,:)

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
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), sinx(step), fftLege(step), symL(19,step), asymL(19,step), &
            & sumLegendre(38,step,0:this%nFourier-1), fft(8,step,0:this%nFourier/2-1), fftC(8,step), cr(4,this%jms1)              )
      
      cr = cmplx(0._dbl, 0._dbl, kind=dbl)

      do i = 1, this%nLegendre, step
        cosx    = this%roots(i:i+step-1)
        sinx    = sqrt(1 - cosx**2)
        fftLege = this%fftLege(i:i+step-1)
        
        pmm = 1._dbl; mj = 0; sumLegendre = cmplx(0._dbl, 0._dbl, kind=dbl)
          do m = 0, this%maxj
            pmj2 = 0._dbl
            pmj1 = 0._dbl
            pmj  = 1._dbl

            symL  = cmplx(0._dbl, 0._dbl, kind=dbl)
            asymL = cmplx(0._dbl, 0._dbl, kind=dbl)

            s = -1

            do j = m, this%maxj
              mj = mj+1 ; s = -s

              if (s == 1) then
                do concurrent ( i2=1:step, i1=1:19 )
                  symL(i1,i2) = symL(i1,i2) + cc(i1,mj) * pmj(i2)
                end do
              else
                do concurrent ( i2=1:step, i1=1:19 )
                  asymL(i1,i2) = asymL(i1,i2) + cc(i1,mj) * pmj(i2)
                end do
              end if

              pmj2 = pmj1
              pmj1 = pmj
              pmj  = this%amjrr(mj+1) * cosx * pmj1 - this%bmjrr(mj+1) * pmj2
            end do

            do concurrent ( i2=1:step, i1=1:19 )
              sumLegendre(i1   ,i2,m) = (symL(i1,i2)+asymL(i1,i2)) * pmm(i2)
              sumLegendre(i1+19,i2,m) = (symL(i1,i2)-asymL(i1,i2)) * pmm(i2)
            end do

            pmm = -this%cmmrr(m) * sinx * pmm; if (maxval(abs(pmm)) < 1.0d-55) exit
          end do

        call fftw_execute_dft(this%fftw_38_forw, sumLegendre, sumLegendre)
        do i1 = 0, this%nFourier/2-1
          do i2 = 1, step
            iL = 2*i1

            fft(1,i2,i1)%re = sumLegendre( 1,i2,iL)%re * sumLegendre( 4,i2,iL)%re + &    
                            & sumLegendre( 2,i2,iL)%re * sumLegendre( 5,i2,iL)%re + &    
                            & sumLegendre( 3,i2,iL)%re * sumLegendre( 6,i2,iL)%re
            fft(2,i2,i1)%re = sumLegendre( 4,i2,iL)%re * sumLegendre( 7,i2,iL)%re + &
                            & sumLegendre( 5,i2,iL)%re * sumLegendre( 8,i2,iL)%re + &
                            & sumLegendre( 6,i2,iL)%re * sumLegendre( 9,i2,iL)%re + &
                            & sumLegendre(16,i2,iL)%re * sumLegendre(19,i2,iL)%re
            fft(3,i2,i1)%re = sumLegendre( 4,i2,iL)%re * sumLegendre(10,i2,iL)%re + &
                            & sumLegendre( 5,i2,iL)%re * sumLegendre(11,i2,iL)%re + &
                            & sumLegendre( 6,i2,iL)%re * sumLegendre(12,i2,iL)%re + &
                            & sumLegendre(17,i2,iL)%re * sumLegendre(19,i2,iL)%re
            fft(4,i2,i1)%re = sumLegendre( 4,i2,iL)%re * sumLegendre(13,i2,iL)%re + &
                            & sumLegendre( 5,i2,iL)%re * sumLegendre(14,i2,iL)%re + &
                            & sumLegendre( 6,i2,iL)%re * sumLegendre(15,i2,iL)%re + &
                            & sumLegendre(18,i2,iL)%re * sumLegendre(19,i2,iL)%re
            fft(5,i2,i1)%re = sumLegendre(20,i2,iL)%re * sumLegendre(23,i2,iL)%re + &    
                            & sumLegendre(21,i2,iL)%re * sumLegendre(24,i2,iL)%re + &    
                            & sumLegendre(22,i2,iL)%re * sumLegendre(25,i2,iL)%re
            fft(6,i2,i1)%re = sumLegendre(23,i2,iL)%re * sumLegendre(26,i2,iL)%re + &
                            & sumLegendre(24,i2,iL)%re * sumLegendre(27,i2,iL)%re + &
                            & sumLegendre(25,i2,iL)%re * sumLegendre(28,i2,iL)%re + &
                            & sumLegendre(35,i2,iL)%re * sumLegendre(38,i2,iL)%re
            fft(7,i2,i1)%re = sumLegendre(23,i2,iL)%re * sumLegendre(29,i2,iL)%re + &
                            & sumLegendre(24,i2,iL)%re * sumLegendre(30,i2,iL)%re + &
                            & sumLegendre(25,i2,iL)%re * sumLegendre(31,i2,iL)%re + &
                            & sumLegendre(36,i2,iL)%re * sumLegendre(38,i2,iL)%re
            fft(8,i2,i1)%re = sumLegendre(23,i2,iL)%re * sumLegendre(32,i2,iL)%re + &
                            & sumLegendre(24,i2,iL)%re * sumLegendre(33,i2,iL)%re + &
                            & sumLegendre(25,i2,iL)%re * sumLegendre(34,i2,iL)%re + &
                            & sumLegendre(37,i2,iL)%re * sumLegendre(38,i2,iL)%re
            
            iL = iL+1
            
            fft(1,i2,i1)%im = sumLegendre( 1,i2,iL)%re * sumLegendre( 4,i2,iL)%re + &    
                            & sumLegendre( 2,i2,iL)%re * sumLegendre( 5,i2,iL)%re + &    
                            & sumLegendre( 3,i2,iL)%re * sumLegendre( 6,i2,iL)%re
            fft(2,i2,i1)%im = sumLegendre( 4,i2,iL)%re * sumLegendre( 7,i2,iL)%re + &
                            & sumLegendre( 5,i2,iL)%re * sumLegendre( 8,i2,iL)%re + &
                            & sumLegendre( 6,i2,iL)%re * sumLegendre( 9,i2,iL)%re + &
                            & sumLegendre(16,i2,iL)%re * sumLegendre(19,i2,iL)%re
            fft(3,i2,i1)%im = sumLegendre( 4,i2,iL)%re * sumLegendre(10,i2,iL)%re + &
                            & sumLegendre( 5,i2,iL)%re * sumLegendre(11,i2,iL)%re + &
                            & sumLegendre( 6,i2,iL)%re * sumLegendre(12,i2,iL)%re + &
                            & sumLegendre(17,i2,iL)%re * sumLegendre(19,i2,iL)%re
            fft(4,i2,i1)%im = sumLegendre( 4,i2,iL)%re * sumLegendre(13,i2,iL)%re + &
                            & sumLegendre( 5,i2,iL)%re * sumLegendre(14,i2,iL)%re + &
                            & sumLegendre( 6,i2,iL)%re * sumLegendre(15,i2,iL)%re + &
                            & sumLegendre(18,i2,iL)%re * sumLegendre(19,i2,iL)%re
            fft(5,i2,i1)%im = sumLegendre(20,i2,iL)%re * sumLegendre(23,i2,iL)%re + &    
                            & sumLegendre(21,i2,iL)%re * sumLegendre(24,i2,iL)%re + &    
                            & sumLegendre(22,i2,iL)%re * sumLegendre(25,i2,iL)%re
            fft(6,i2,i1)%im = sumLegendre(23,i2,iL)%re * sumLegendre(26,i2,iL)%re + &
                            & sumLegendre(24,i2,iL)%re * sumLegendre(27,i2,iL)%re + &
                            & sumLegendre(25,i2,iL)%re * sumLegendre(28,i2,iL)%re + &
                            & sumLegendre(35,i2,iL)%re * sumLegendre(38,i2,iL)%re
            fft(7,i2,i1)%im = sumLegendre(23,i2,iL)%re * sumLegendre(29,i2,iL)%re + &
                            & sumLegendre(24,i2,iL)%re * sumLegendre(30,i2,iL)%re + &
                            & sumLegendre(25,i2,iL)%re * sumLegendre(31,i2,iL)%re + &
                            & sumLegendre(36,i2,iL)%re * sumLegendre(38,i2,iL)%re
            fft(8,i2,i1)%im = sumLegendre(23,i2,iL)%re * sumLegendre(32,i2,iL)%re + &
                            & sumLegendre(24,i2,iL)%re * sumLegendre(33,i2,iL)%re + &
                            & sumLegendre(25,i2,iL)%re * sumLegendre(34,i2,iL)%re + &
                            & sumLegendre(37,i2,iL)%re * sumLegendre(38,i2,iL)%re
          end do
        end do
        call fftw_execute_dft(this%fftw_08_back, fft, fft); fft(:,:,0)%im = 0._dbl

        pmm = 1._dbl; mj = 1; mult = exp(-2 * pi * cunit / this%nFourier)
          do m = 0, this%jmax+1 
            pmj2 = 0._dbl
            pmj1 = 0._dbl
            pmj  = 1._dbl

            if (m == 0) then
              jexp = cunit
                fftC = ( (1-cunit) * fft(:,:,0) + (1+cunit) * conjg(fft(:,:,0)) ) / 2
            else
              jexp = mult * jexp
                fftC = ( (1-jexp) * fft(:,:,m) + (1+jexp) * conjg(fft(:,:,this%nFourier/2-m)) ) / 2
            end if

            do concurrent (i2=1:step, i1=1:8)
              fftC(i1,i2) = fftLege(i2) * fftC(i1,i2) * pmm(i2)
            end do

            do j = m, this%jmax+1
              do i2 = 1, step
                do i1 = 1, 4
                  cr(i1,mj) = cr(i1,mj) + pmj(i2) * ( fftC(i1,i2) + (-1)**(j+m) * fftC(i1+4,i2) )
                end do
              end do

              mj = mj+1
                pmj2 = pmj1
                pmj1 = pmj
                pmj  = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
            end do
            
            pmm = -this%cmmrr(m) * sinx * pmm; if (maxval(abs(pmm)) < 1.0d-55) exit
          end do
      end do

      cr = cr / 4 / this%nLegendre**2 / this%nFourier / sqrt(pi)

    deallocate( cc, fft, fftC, sumLegendre, pmm, pmj, pmj1, pmj2, cosx, sinx, fftLege, symL, asymL )
    
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

    call fftw_destroy_plan( this%fftw_38_forw )
    call fftw_destroy_plan( this%fftw_08_back )

  end subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub

end submodule vcsv_vcvv_vcvgv