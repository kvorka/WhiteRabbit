submodule (SphericalHarmonics) vcsv_vcvgv
  implicit none

  contains

  subroutine init_vcsv_vcvgv_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer,              pointer       :: ip(:) => null()
    real(kind=dbl),       allocatable   :: field_re(:,:,:)
    complex(kind=dbl),    allocatable   :: field(:,:,:)

    allocate(field(16,step,this%nFourier/2+1), field_re(16,step,this%nFourier))
      this%fftw_16_c2r = fftw_plan_many_dft_c2r( 1, [this%nFourier], 16*step, field,    ip, 16*step, 1,            &
                                               &                              field_re, ip, 16*step, 1, fftw_flags )
    deallocate(field, field_re)

    allocate(field_re(3,step,this%nFourier), field(step,3,this%nFourier/2+1))
      this%fftw_03_r2c = fftw_plan_many_dft_r2c( 1, [this%nFourier], 3*step, field_re, ip, 3*step, 1,            &
                                                &                            field   , ip, 3*step, 1, fftw_flags )
    deallocate(field, field_re)

    write(*,*) 'vcsv_vcvgv initialized'

  end subroutine init_vcsv_vcvgv_sub

  function vcsv_vcvgv_fn(this, ri, dv_r, v) result(cjml)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(:), v(:)
    complex(kind=dbl)                 :: cjml(this%jmv)
    integer                           :: i, j, m, l, ijm, lm, ijml, i1, i2, mj, s, iL, lm1, lm2
    real(kind=dbl)                    :: fac
    real(kind=dbl),       allocatable :: gc(:), pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), fftLege(:), sinx(:)
    real(kind=dbl),       allocatable :: grid(:,:,:), fft(:,:,:)
    complex(kind=dbl)                 :: cr12
    complex(kind=dbl),    allocatable :: sum1(:), sum2(:), sum3(:), cab(:,:), cc(:,:), cr(:,:)
    complex(kind=dbl),    allocatable :: sumLegendreN(:,:,:), sumLegendreS(:,:,:), fftNC(:,:,:), fftSC(:,:,:)

    allocate(cab(5,this%jmv1), sum1(2), sum2(2), sum3(2), gc(2)) ; cab = czero

      do ijml = 1, this%jmv
        cab(1,ijml) = v(ijml)
        cab(5,ijml) = dv_r(ijml)
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
            cab(2,ijml) =         ( sum1(1) - sum2(1) ) / ri
            cab(3,ijml) = cunit * ( sum1(1) + sum2(1) ) / ri
            cab(4,ijml) =         ( sum3(1)           ) / ri

          ijml = 3*(j*(j+1)/2+m)+(j+1)-j
            cab(2,ijml) =         ( sum1(2) - sum2(2) ) / ri
            cab(3,ijml) = cunit * ( sum1(2) + sum2(2) ) / ri
            cab(4,ijml) =         ( sum3(2)           ) / ri
        end do
      end do

    deallocate(sum1, sum2, sum3, gc)
    allocate(cc(16,this%jms2), sum1(5), sum2(5), sum3(5)) ; cc = czero
    
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
            do i1 = 1, 5
              cc(3*(i1-1)+1, mj) =           sum1(i1) - sum2(i1)
              cc(3*(i1-1)+2, mj) = cunit * ( sum1(i1) + sum2(i1) )
              cc(3*(i1-1)+3, mj) =           sum3(i1)
            end do

            if (j <= this%jmax) then
              ijml = 3*(j*(j+1)/2+m)

              if (j == 0) then
                cc(16,mj) = -v(ijml+1)
              else
                cc(16,mj) = sqrt(j/(2*j+1._dbl)) * v(ijml-1) - sqrt((j+1)/(2*j+1._dbl)) * v(ijml+1)
              end if
            end if
        end do
      end do
      
      do mj = 1, this%jms2
        do i1 = 1, 2
          cc(i1,mj) = cc(i1,mj) / 2
        end do

        do i1 = 4, 15
          cc(i1,mj) = cc(i1,mj) / 2
        end do
      end do

    deallocate(cab, sum1, sum2, sum3)
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), fftLege(step), sinx(step), &
            & sumLegendreN(16,step,0:this%nFourier/2), sumLegendreS(16,step,0:this%nFourier/2),    &                                 
            & cr(3,this%jms1), grid(16,step,0:this%nFourier-1), fft(3,step,0:this%nFourier-1),     & 
            & fftNC(3,step,0:this%nFourier/2), fftSC(3,step,0:this%nFourier/2)                     )
      
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
          
          j = m
            s = +1 ; mj = mj+1 
            
            do concurrent (i2=1:step, i1=1:16)
              sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj)
              sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj)
            end do
            
          do j = m+1, this%maxj
            s = -s ; mj = mj+1

            pmj2 = pmj1
            pmj1 = pmj
            pmj  = ( cosx * pmj1 - this%ish(mj-1) * pmj2) / this%ish(mj)

            do concurrent (i2=1:step, i1=1:16)
              sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj) * pmj(i2)
              sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj) * pmj(i2) * s
            end do
          end do
          
          do concurrent ( i2=1:step, i1=1:16 )
            sumLegendreN(i1,i2,0)%im = 0._dbl
            sumLegendreS(i1,i2,0)%im = 0._dbl
          end do
          
        do m = 1, this%maxj
          fac = sqrt(( 2*m+1._dbl ) / ( 2*m ))
          pmm = -fac * sinx * pmm
          if (maxval(abs(pmm)) < 1.0d-55) exit
          
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl

          j = m
          s = +1 ; mj = mj+1
            
          do concurrent (i2=1:step, i1=1:16)
            sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj)
            sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj)
          end do
          
          do j = m+1, this%maxj
            s = -s ; mj = mj+1

            pmj2 = pmj1
            pmj1 = pmj
            pmj  = ( cosx * pmj1 - this%ish(mj-1) * pmj2) / this%ish(mj)
            
            do concurrent (i2=1:step, i1=1:16)
              sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj) * pmj(i2)
              sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj) * pmj(i2) * s
            end do
          end do
          
          do concurrent ( i2=1:step, i1=1:16 )
            sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) * pmm(i2)
            sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) * pmm(i2)
          end do
        end do

        call fftw_execute_dft_c2r(this%fftw_16_c2r, sumLegendreN, grid)
        
        do concurrent (i1=0:this%nFourier-1, i2=1:step)
          fft(1,i2,i1) = grid( 1,i2,i1) * grid( 4,i2,i1) + &
                       & grid( 2,i2,i1) * grid( 5,i2,i1) + &
                       & grid( 3,i2,i1) * grid( 6,i2,i1) + &
                       & grid(13,i2,i1) * grid(16,i2,i1)
          fft(2,i2,i1) = grid( 1,i2,i1) * grid( 7,i2,i1) + &
                       & grid( 2,i2,i1) * grid( 8,i2,i1) + &
                       & grid( 3,i2,i1) * grid( 9,i2,i1) + &
                       & grid(14,i2,i1) * grid(16,i2,i1)
          fft(3,i2,i1) = grid( 1,i2,i1) * grid(10,i2,i1) + &
                       & grid( 2,i2,i1) * grid(11,i2,i1) + &
                       & grid( 3,i2,i1) * grid(12,i2,i1) + &
                       & grid(15,i2,i1) * grid(16,i2,i1)
        end do

        call fftw_execute_dft_r2c(this%fftw_03_r2c, fft, fftNC)

        call fftw_execute_dft_c2r(this%fftw_16_c2r, sumLegendreS, grid)
        
        do concurrent (i1=0:this%nFourier-1, i2=1:step)
          fft(1,i2,i1) = grid( 1,i2,i1) * grid( 4,i2,i1) + &
                       & grid( 2,i2,i1) * grid( 5,i2,i1) + &
                       & grid( 3,i2,i1) * grid( 6,i2,i1) + &
                       & grid(13,i2,i1) * grid(16,i2,i1)
          fft(2,i2,i1) = grid( 1,i2,i1) * grid( 7,i2,i1) + &
                       & grid( 2,i2,i1) * grid( 8,i2,i1) + &
                       & grid( 3,i2,i1) * grid( 9,i2,i1) + &
                       & grid(14,i2,i1) * grid(16,i2,i1)
          fft(3,i2,i1) = grid( 1,i2,i1) * grid(10,i2,i1) + &
                       & grid( 2,i2,i1) * grid(11,i2,i1) + &
                       & grid( 3,i2,i1) * grid(12,i2,i1) + &
                       & grid(15,i2,i1) * grid(16,i2,i1)
        end do

        call fftw_execute_dft_r2c(this%fftw_03_r2c, fft, fftSC)
        
        pmm  = 1._dbl
        mj   = 0
        
        m = 0
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl
          
          do concurrent (i2=1:step, i1=1:3)
            fftNC(i1,i2,m) = fftLege(i2) * fftNC(i1,i2,m) ; fftNC(i1,i2,m)%im = 0._dbl
            fftSC(i1,i2,m) = fftLege(i2) * fftSC(i1,i2,m) ; fftSC(i1,i2,m)%im = 0._dbl
          end do
          
          j = m
            s = +1 ; mj = mj+1
            
            do concurrent(i2=1:step, i1=1:3)
              cr(i1,mj) = cr(i1,mj) +  fftNC(i1,i2,m) + fftSC(i1,i2,m)
            end do
          
          do j = m+1, this%jmax+1
            s = -s ; mj = mj+1
            
            pmj2 = pmj1
            pmj1 = pmj
            pmj  = ( cosx * pmj1 - this%ish(mj+m-1) * pmj2) / this%ish(mj+m)
            
            do concurrent(i2=1:step, i1=1:3)
              cr(i1,mj) = cr(i1,mj) + pmj(i2) * ( fftNC(i1,i2,m) + s * fftSC(i1,i2,m) )
            end do
          end do
        
        do m = 1, this%jmax+1
          fac = sqrt(( 2*m+1._dbl ) / ( 2*m ))
          pmm = -fac * sinx * pmm
          if (maxval(abs(pmm)) < 1.0d-55) exit
          
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl

          do concurrent (i2=1:step, i1=1:3)
            fftNC(i1,i2,m) = fftLege(i2) * fftNC(i1,i2,m) * pmm(i2)
            fftSC(i1,i2,m) = fftLege(i2) * fftSC(i1,i2,m) * pmm(i2)
          end do
          
          j = m
            s = +1 ; mj = mj+1
            
            do concurrent(i2=1:step, i1=1:3)
              cr(i1,mj) = cr(i1,mj) + fftNC(i1,i2,m) + fftSC(i1,i2,m)
            end do
          
          do j = m+1, this%jmax+1
            s = -s ; mj = mj+1
            
            pmj2 = pmj1
            pmj1 = pmj
            pmj  = ( cosx * pmj1 - this%ish(mj+m-1) * pmj2) / this%ish(mj+m)
            
            do concurrent(i2=1:step, i1=1:3)
              cr(i1,mj) = cr(i1,mj) + pmj(i2) * ( fftNC(i1,i2,m) + s * fftSC(i1,i2,m) )
            end do
          end do
        end do
      end do

    deallocate( cc, sumLegendreN, sumLegendreS, grid, fft, fftNC, fftSC, pmm, pmj, pmj1, pmj2, cosx, sinx, fftLege )
    
    fac = 1 / ( 4 * this%nLegendre**2 * this%nFourier * sqrt(pi) )
    
    do mj = 1, this%jms1
      cr(1,mj) = cr(1,mj) * fac / 2
      cr(2,mj) = cr(2,mj) * fac / 2 * cunit
      cr(3,mj) = cr(3,mj) * fac

      cr12     = cr(1,mj) - cr(2,mj)
      cr(2,mj) = cr(1,mj) + cr(2,mj)
      cr(1,mj) = cr12
    end do
    
    j = 0
      m = 0
        ijml = -1
        lm   = m*(this%maxj)-m*(m+1)/2+j+1
        lm2  = lm+this%maxj-m
        
        cjml(ijml+2) =        cr(1,lm2 )   * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                     &        cr(3,lm+1)   * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                     & conjg( cr(1,lm2 ) ) * cleb1_fn(j+1,m-1,1,+1,j,m) ; cjml(ijml+2)%im = 0._dbl

    do j = 1, this%jmax
      m = 0
        ijml = ijml+3
        lm   = m*(this%maxj)-m*(m+1)/2+j+1
        lm2  = lm+this%maxj-m-1
        
        cjml(ijml  ) =        cr(1,lm2-1)  * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                     &        cr(3,lm -1)  * cleb1_fn(j-1,m+0,1, 0,j,m) + &
                     & conjg( cr(1,lm2-1) )* cleb1_fn(j-1,m-1,1,+1,j,m)
        cjml(ijml+1) =        cr(1,lm2  )  * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                     &        cr(3,lm   )  * cleb1_fn(j  ,m+0,1, 0,j,m) + &
                     & conjg( cr(1,lm2  ) )* cleb1_fn(j  ,m-1,1,+1,j,m)
        cjml(ijml+2) =        cr(1,lm2+1)  * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                     &        cr(3,lm +1)  * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                     & conjg( cr(1,lm2+1) )* cleb1_fn(j+1,m-1,1,+1,j,m)
        
        cjml(ijml  )%im = 0._dbl
        cjml(ijml+1)%re = 0._dbl
        cjml(ijml+2)%im = 0._dbl

      do m = 1, j
        ijml = ijml+3
        lm   = m*(this%maxj)-m*(m+1)/2+j+1
        lm1  = lm - this%maxj + m
        lm2  = lm + this%maxj - m - 1
        
        cjml(ijml  ) = cr(1,lm2-1) * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                     & cr(3,lm -1) * cleb1_fn(j-1,m+0,1, 0,j,m) - &
                     & cr(2,lm1-1) * cleb1_fn(j-1,m-1,1,+1,j,m)
        cjml(ijml+1) = cr(1,lm2  ) * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                     & cr(3,lm   ) * cleb1_fn(j  ,m+0,1, 0,j,m) - &
                     & cr(2,lm1  ) * cleb1_fn(j  ,m-1,1,+1,j,m)
        cjml(ijml+2) = cr(1,lm2+1) * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                     & cr(3,lm +1) * cleb1_fn(j+1,m+0,1, 0,j,m) - &
                     & cr(2,lm1+1) * cleb1_fn(j+1,m-1,1,+1,j,m)
      end do
    end do

    deallocate(cr)

  end function vcsv_vcvgv_fn

  subroutine deallocate_fftw_vcsv_vcvgv_sub(this)
    class(T_lateralGrid), intent(inout) :: this

    call fftw_destroy_plan( this%fftw_16_c2r )
    call fftw_destroy_plan( this%fftw_03_r2c )

  end subroutine deallocate_fftw_vcsv_vcvgv_sub

end submodule vcsv_vcvgv