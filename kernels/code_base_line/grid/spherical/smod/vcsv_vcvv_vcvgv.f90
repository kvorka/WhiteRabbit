submodule (SphericalHarmonics) vcsv_vcvv_vcvgv
  implicit none

  contains

  subroutine init_vcsv_vcvv_vcvgv_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer,              pointer       :: ip(:) => null()
    real(kind=dbl),       allocatable   :: field_re(:,:,:)
    complex(kind=dbl),    allocatable   :: field(:,:,:)

    allocate(field(19,step,this%nFourier/2+1), field_re(19,step,this%nFourier))
      this%fftw_19_c2r = fftw_plan_many_dft_c2r( 1, [this%nFourier], 19*step, field,    ip, 19*step, 1,            &
                                               &                              field_re, ip, 19*step, 1, fftw_flags )
    deallocate(field, field_re)

    allocate(field(4,step,this%nFourier/2+1), field_re(4,step,this%nFourier))
      this%fftw_04_r2c = fftw_plan_many_dft_r2c( 1, [this%nFourier], 4*step, field_re, ip, 4*step, 1,            &
                                               &                             field   , ip, 4*step, 1, fftw_flags )
    deallocate(field, field_re)

    write(*,*) 'vcsv_vcvv_vcvgv initialized'

  end subroutine init_vcsv_vcvv_vcvgv_sub

  subroutine vcsv_vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(:), q(:), v(:)
    complex(kind=dbl),    intent(out) :: cjm(:,:)
    integer                           :: i, j, m, l, ijm, lm, ijml, i1, i2, mj, s, iL, lm1, lm2, jml_int
    real(kind=dbl)                    :: fac
    real(kind=dbl),       allocatable :: gc(:), pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), fftLege(:)
    real(kind=dbl),       allocatable :: grid(:,:,:), fft(:,:,:)
    complex(kind=dbl),    allocatable :: sum1(:), sum2(:), sum3(:)
    complex(kind=dbl),    allocatable :: cab(:,:), cc(:,:), cr(:,:), symL(:,:), asymL(:,:)
    complex(kind=dbl),    allocatable :: sumLegendreN(:,:,:), sumLegendreS(:,:,:), fftNC(:,:,:), fftSC(:,:,:)

    allocate(cab(6,this%jmv1), sum1(2), sum2(2), sum3(2), gc(2)) ; cab = czero

      do ijml = 1, this%jmv
        cab(1,ijml) = q(ijml)
        cab(2,ijml) = v(ijml)
        cab(6,ijml) = dv_r(ijml)
      end do

      do j = 1, this%jmax+1
        gc(1) = sqrt(j*(j+1._dbl)*(j+1._dbl)/(2*j+1._dbl))
        gc(2) = sqrt(j*(j       )*(j+1._dbl)/(2*j+1._dbl))

        do m = 0, j
          sum1 = czero
          sum2 = czero
          sum3 = czero

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
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), fftLege(step), sinx(step),                             &
            & symL(19,step), asymL(19,step), sumLegendreN(19,step,0:this%nFourier/2), sumLegendreS(19,step,0:this%nFourier/2), &                                 
            & cr(4,this%jms1), grid(19,step,0:this%nFourier-1), fft(4,step,0:this%nFourier-1),                                 & 
            & fftNC(4,step,0:this%nFourier/2), fftSC(4,step,0:this%nFourier/2))
      
      cr = czero

      do i = 1, this%nLegendre, step
        cosx    = this%roots(i:i+step-1)
        sinx    = sqrt(1-cosx**2)
        fftLege = this%fftLege(i:i+step-1)
        
        pmm = 1._dbl
        mj  = 0

        sumLegendreN = czero
        sumLegendreS = czero

        m = 0
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl
          
          symL  = czero
          asymL = czero
          
          j = m
            mj = mj+1 
            
            do concurrent ( i2=1:step, i1=1:19 )
              symL(i1,i2) = symL(i1,i2) + cc(i1,mj)
            end do
            
          do j = 1, (this%maxj-m)/2
            mj = mj+2
            
            pmj2 = pmj1
            pmj1 = pmj
            pmj = ( cosx * pmj1 - this%ish(mj-2) * pmj2 ) / this%ish(mj-1)

            do concurrent ( i2=1:step, i1=1:19 )
              asymL(i1,i2) = asymL(i1,i2) + cc(i1,mj-1) * pmj(i2)
            end do

            pmj2 = pmj1
            pmj1 = pmj
            pmj = ( cosx * pmj1 - this%ish(mj-1) * pmj2 ) / this%ish(mj)

            do concurrent ( i2=1:step, i1=1:19 )
              symL(i1,i2) = symL(i1,i2) + cc(i1,mj) * pmj(i2)
            end do
          end do

          if (mod((this%maxj-m),2) /= 0) then
            mj = mj+1

            pmj2 = pmj1
            pmj1 = pmj
            pmj = ( cosx * pmj1 - this%ish(mj-1) * pmj2 ) / this%ish(mj)

            do concurrent ( i2=1:step, i1=1:19 )
              asymL(i1,i2) = asymL(i1,i2) + cc(i1,mj) * pmj(i2)
            end do
          end if
          
          do concurrent ( i2=1:step, i1=1:19 )
            sumLegendreN(i1,i2,0)%re = symL(i1,i2)%re+asymL(i1,i2)%re ; sumLegendreN(i1,i2,0)%im = 0._dbl
            sumLegendreS(i1,i2,0)%re = symL(i1,i2)%re-asymL(i1,i2)%re ; sumLegendreS(i1,i2,0)%im = 0._dbl
          end do
          
        do m = 1, this%maxj
          fac = -sqrt( ( 2*m+1 ) / ( 2._dbl*m ) )
          pmm = fac * sinx * pmm
          if (maxval(abs(pmm)) < 1.0d-55) exit
          
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl
          
          symL  = czero
          asymL = czero

          j = m
            mj = mj+1
            
            do concurrent ( i2=1:step, i1=1:19 )
              symL(i1,i2) = symL(i1,i2) + cc(i1,mj)
            end do
          
          do j = 1, (this%maxj-m)/2
            mj = mj+2

            pmj2 = pmj1
            pmj1 = pmj
            pmj = ( cosx * pmj1 - this%ish(mj-2) * pmj2) / this%ish(mj-1)

            do concurrent ( i2=1:step, i1=1:19 )
              asymL(i1,i2) = asymL(i1,i2) + cc(i1,mj-1) * pmj(i2)
            end do

            pmj2 = pmj1
            pmj1 = pmj
            pmj = ( cosx * pmj1 - this%ish(mj-1) * pmj2) / this%ish(mj)

            do concurrent ( i2=1:step, i1=1:19 )
              symL(i1,i2) = symL(i1,i2) + cc(i1,mj) * pmj(i2)
            end do
          end do

          if (mod((this%maxj-m),2) /= 0) then
            mj = mj+1

            pmj2 = pmj1
            pmj1 = pmj
            pmj = ( cosx * pmj1 - this%ish(mj-1) * pmj2) / this%ish(mj)

            do concurrent ( i2=1:step, i1=1:19 )
              asymL(i1,i2) = asymL(i1,i2) + cc(i1,mj) * pmj(i2)
            end do
          end if
          
          do concurrent ( i2=1:step, i1=1:19 )
            sumLegendreN(i1,i2,m) = (symL(i1,i2)+asymL(i1,i2)) * pmm(i2)
            sumLegendreS(i1,i2,m) = (symL(i1,i2)-asymL(i1,i2)) * pmm(i2)
          end do
        end do

        call fftw_execute_dft_c2r(this%fftw_19_c2r, sumLegendreN, grid)
        
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

        call fftw_execute_dft_r2c(this%fftw_04_r2c, fft, fftNC)

        call fftw_execute_dft_c2r(this%fftw_19_c2r, sumLegendreS, grid)
        
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
        
        call fftw_execute_dft_r2c(this%fftw_04_r2c, fft, fftSC)
        
        pmm  = 1._dbl
        mj   = 0
        
        m = 0
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl
          
          do concurrent (i2=1:step, i1=1:4)
            fftNC(i1,i2,m) = fftLege(i2) * fftNC(i1,i2,m) ; fftNC(i1,i2,m)%im = 0._dbl
            fftSC(i1,i2,m) = fftLege(i2) * fftSC(i1,i2,m) ; fftSC(i1,i2,m)%im = 0._dbl
          end do
          
          j = m
            s = +1 ; mj = mj+1

            do concurrent (i2=1:step, i1=1:4)
              cr(i1,mj) = cr(i1,mj) + fftNC(i1,i2,m) + fftSC(i1,i2,m)
            end do
          
          do j = m+1, this%jmax+1
            s = -s ; mj = mj+1
            
            pmj2 = pmj1
            pmj1 = pmj
            pmj = ( cosx * pmj1 - this%ish(mj+m-1) * pmj2) / this%ish(mj+m)
            
            do concurrent (i2=1:step, i1=1:4)
              cr(i1,mj) = cr(i1,mj) + pmj(i2) * ( fftNC(i1,i2,m) + s * fftSC(i1,i2,m) )
            end do
          end do
        
        do m = 1, this%jmax+1
          fac = -sqrt( ( 2*m+1 ) / ( 2._dbl*m ) )
          pmm = fac * sinx * pmm
          if (maxval(abs(pmm)) < 1.0d-55) exit
          
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl

            do concurrent (i2=1:step, i1=1:4)
              fftNC(i1,i2,m) = fftLege(i2) * fftNC(i1,i2,m) * pmm(i2)
              fftSC(i1,i2,m) = fftLege(i2) * fftSC(i1,i2,m) * pmm(i2)
            end do
          
          j = m
            s = +1 ; mj = mj+1
            
            do concurrent (i2=1:step, i1=1:4)
              cr(i1,mj) = cr(i1,mj) + fftNC(i1,i2,m) + fftSC(i1,i2,m)
            end do
          
          do j = m+1, this%jmax+1
            s = -s ; mj = mj+1
            
            pmj2 = pmj1
            pmj1 = pmj
            pmj = ( cosx * pmj1 - this%ish(mj+m-1) * pmj2) / this%ish(mj+m)
            
            do concurrent (i2=1:step, i1=1:4)
              cr(i1,mj) = cr(i1,mj) + pmj(i2) * ( fftNC(i1,i2,m) + s * fftSC(i1,i2,m) )
            end do
          end do
        end do
      end do

    deallocate( cc, sumLegendreN, sumLegendreS, grid, fft, fftNC, fftSC, pmm, pmj, pmj1, pmj2, cosx, sinx, fftLege, symL, asymL )
    
    fac = 1 / ( 4 * this%nLegendre**2 * this%nFourier * sqrt(pi) )
    
    do mj = 1, this%jms1 
      cr(1,mj) = cr(1,mj) * fac
      cr(2,mj) = cr(2,mj) * fac / 2
      cr(3,mj) = cr(3,mj) * fac / 2 * cunit
      cr(4,mj) = cr(4,mj) * fac
    end do
    
    j = 0
      m = 0
        ijm = 1
        lm  = m*(this%maxj)-m*(m+1)/2+j+1
        lm2 = lm+this%maxj-m

        cjm(1,ijm) = cr(1,lm)
        cjm(2,ijm) = czero
        cjm(3,ijm) = czero
        cjm(4,ijm) = (        cr(2,lm2 ) - cr(3,lm2)   ) * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                   & (        cr(4,lm+1)               ) * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                   & ( conjg( cr(2,lm2 ) - cr(3,lm2) ) ) * cleb1_fn(j+1,m-1,1,+1,j,m)
        
        cjm(1,ijm)%im = 0._dbl ; cjm(4,ijm)%im = 0._dbl

    do j = 1, this%jmax
      m = 0
        ijm = ijm+1
        lm  = m*(this%maxj)-m*(m+1)/2+j+1
        lm2 = lm+this%maxj-m-1

        cjm(1,ijm) = cr(1,lm)
        cjm(2,ijm) = (        cr(2,lm2-1) - cr(3,lm2-1)   ) * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                   & (        cr(4,lm -1)                 ) * cleb1_fn(j-1,m+0,1, 0,j,m) + &
                   & ( conjg( cr(2,lm2-1) - cr(3,lm2-1) ) ) * cleb1_fn(j-1,m-1,1,+1,j,m)
        cjm(3,ijm) = (        cr(2,lm2  ) - cr(3,lm2  )   ) * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                   & (        cr(4,lm   )                 ) * cleb1_fn(j  ,m+0,1, 0,j,m) + &
                   & ( conjg( cr(2,lm2  ) - cr(3,lm2  ) ) ) * cleb1_fn(j  ,m-1,1,+1,j,m)
        cjm(4,ijm) = (        cr(2,lm2+1) - cr(3,lm2+1)   ) * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                   & (        cr(4,lm +1)                 ) * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                   & ( conjg( cr(2,lm2+1) - cr(3,lm2+1) ) ) * cleb1_fn(j+1,m-1,1,+1,j,m)
        
        cjm(1,ijm)%im = 0._dbl ; cjm(2,ijm)%im = 0._dbl ; cjm(3,ijm)%re = 0._dbl ; cjm(4,ijm)%im = 0._dbl

      do m = 1, j
        ijm = ijm+1
        lm  = m*(this%maxj)-m*(m+1)/2+j+1
        lm1 = lm - this%maxj + m
        lm2 = lm + this%maxj - m - 1

        cjm(1,ijm) = cr(1,lm)
        cjm(2,ijm) = ( cr(2,lm2-1) - cr(3,lm2-1) ) * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                   & ( cr(4,lm -1)               ) * cleb1_fn(j-1,m+0,1, 0,j,m) - &
                   & ( cr(2,lm1-1) + cr(3,lm1-1) ) * cleb1_fn(j-1,m-1,1,+1,j,m)
        cjm(3,ijm) = ( cr(2,lm2  ) - cr(3,lm2  ) ) * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                   & ( cr(4,lm   )               ) * cleb1_fn(j  ,m+0,1, 0,j,m) - &
                   & ( cr(2,lm1  ) + cr(3,lm1  ) ) * cleb1_fn(j  ,m-1,1,+1,j,m)
        cjm(4,ijm) = ( cr(2,lm2+1) - cr(3,lm2+1) ) * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                   & ( cr(4,lm +1)               ) * cleb1_fn(j+1,m+0,1, 0,j,m) - &
                   & ( cr(2,lm1+1) + cr(3,lm1+1) ) * cleb1_fn(j+1,m-1,1,+1,j,m)
      end do
    end do

    deallocate(cr)
    
  end subroutine vcsv_vcvv_vcvgv_sub

  subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub(this)
    class(T_lateralGrid), intent(inout) :: this

    call fftw_destroy_plan( this%fftw_19_c2r )
    call fftw_destroy_plan( this%fftw_04_r2c )

  end subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub

end submodule vcsv_vcvv_vcvgv