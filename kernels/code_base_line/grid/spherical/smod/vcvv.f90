submodule (SphericalHarmonics) vcvv
  implicit none

  contains

  subroutine init_vcvv_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer,              pointer       :: ip(:) => null()
    real(kind=dbl),       allocatable   :: field_re(:,:,:), field1_re(:,:)
    complex(kind=dbl),    allocatable   :: field(:,:,:), field1(:,:)

    allocate(field(6,step,this%nFourier/2+1), field_re(6,step,this%nFourier))
      this%fftw_06_c2r = fftw_plan_many_dft_c2r( 1, [this%nFourier], 6*step, field,    ip, 6*step, 1,            &
                                               &                             field_re, ip, 6*step, 1, fftw_flags )
    deallocate(field, field_re)

    allocate(field1_re(step,this%nFourier), field1(step,this%nFourier/2+1))
      this%fftw_01_r2c = fftw_plan_many_dft_r2c( 1, [this%nFourier], step, field1_re, ip, step, 1,            &
                                               &                           field1   , ip, step, 1, fftw_flags )
    deallocate(field1_re, field1)

    write(*,*) 'vcvv initialized'

  end subroutine init_vcvv_sub

  function vcvv_fn(this, cajml, cbjml) result(cjm)
    class(T_lateralGrid), intent(in) :: this
    complex(kind=dbl),    intent(in) :: cajml(:), cbjml(:)
    complex(kind=dbl)                :: cjm(this%jms)
    integer                          :: i, k, j, m, l, mj, i1, i2, s
    real(kind=dbl)                   :: cleb1, cleb2, cleb3, fac
    real(kind=dbl),     allocatable  :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), fftLege(:), fft(:,:), grid(:,:,:), sinx(:)
    complex(kind=dbl),  allocatable  :: cc(:,:), fftNC(:,:), fftSC(:,:), sumLegendreN(:,:,:), sumLegendreS(:,:,:), cr(:)

    cjm = czero
    
    allocate( cc(6,this%jms2) ); cc = czero
    
    mj = 0
      do m = 0, this%maxj
        do j = m, this%maxj
          mj = mj + 1

          do l = abs(j-1), min(this%jmax, j+1)
            if (m == 0) then
              cleb1 = cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(l+j)
              cleb2 = cleb1_fn(j,m,1,+1,l,m+1)
              cleb3 = cleb1_fn(j,m,1, 0,l,m  )

              cc(1,mj) = cc(1,mj) + conjg(cajml(3*(l*(l+1)/2+m+1)+j-l)) * cleb1
              cc(2,mj) = cc(2,mj) +       cajml(3*(l*(l+1)/2+m+1)+j-l)  * cleb2
              cc(3,mj) = cc(3,mj) +       cajml(3*(l*(l+1)/2+m  )+j-l)  * cleb3
              cc(4,mj) = cc(4,mj) + conjg(cbjml(3*(l*(l+1)/2+m+1)+j-l)) * cleb1
              cc(5,mj) = cc(5,mj) +       cbjml(3*(l*(l+1)/2+m+1)+j-l)  * cleb2
              cc(6,mj) = cc(6,mj) +       cbjml(3*(l*(l+1)/2+m  )+j-l)  * cleb3

            else if (l > m+0) then
              cleb1 = cleb1_fn(j,m,1,-1,l,m-1)
              cleb2 = cleb1_fn(j,m,1,+1,l,m+1)
              cleb3 = cleb1_fn(j,m,1, 0,l,m  )

              cc(1,mj) = cc(1,mj) + cajml(3*(l*(l+1)/2+m-1)+j-l) * cleb1
              cc(2,mj) = cc(2,mj) + cajml(3*(l*(l+1)/2+m+1)+j-l) * cleb2
              cc(3,mj) = cc(3,mj) + cajml(3*(l*(l+1)/2+m  )+j-l) * cleb3
              cc(4,mj) = cc(4,mj) + cbjml(3*(l*(l+1)/2+m-1)+j-l) * cleb1
              cc(5,mj) = cc(5,mj) + cbjml(3*(l*(l+1)/2+m+1)+j-l) * cleb2
              cc(6,mj) = cc(6,mj) + cbjml(3*(l*(l+1)/2+m  )+j-l) * cleb3

            else if (l > m-1) then
              cleb1 = cleb1_fn(j,m,1,-1,l,m-1)
              cleb3 = cleb1_fn(j,m,1, 0,l,m  )

              cc(1,mj) = cc(1,mj) + cajml(3*(l*(l+1)/2+m-1)+j-l) * cleb1
              cc(3,mj) = cc(3,mj) + cajml(3*(l*(l+1)/2+m  )+j-l) * cleb3
              cc(4,mj) = cc(4,mj) + cbjml(3*(l*(l+1)/2+m-1)+j-l) * cleb1
              cc(6,mj) = cc(6,mj) + cbjml(3*(l*(l+1)/2+m  )+j-l) * cleb3
          
            else
              cleb1 = cleb1_fn(j,m,1,-1,l,m-1)

              cc(1,mj) = cc(1,mj) + cajml(3*(l*(l+1)/2+m-1)+j-l) * cleb1
              cc(4,mj) = cc(4,mj) + cbjml(3*(l*(l+1)/2+m-1)+j-l) * cleb1
            end if
          end do

          cc(1,mj) =           cc(1,mj)-cc(2,mj)          
          cc(2,mj) = cunit * ( cc(1,mj)+cc(2,mj)+cc(2,mj) )
          cc(3,mj) =     2 *   cc(3,mj)               
          cc(4,mj) =           cc(4,mj)-cc(5,mj)          
          cc(5,mj) = cunit * ( cc(4,mj)+cc(5,mj)+cc(5,mj) )
        end do
      end do
    
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), fftLege(step), sumLegendreN(6,step,0:this%nFourier/2), &
            & sumLegendreS(6,step,0:this%nFourier/2), fft(step,0:this%nFourier-1), grid(6,step,0:this%nFourier-1),             &
            & fftNC(step,0:this%nFourier/2), fftSC(step,0:this%nFourier/2), sinx(step), cr(this%jms1)                          )
            
    cr = czero

    do i = 1, this%nLegendre, step
      cosx    = this%roots(i:i+step-1)
      sinx    = sqrt(1-cosx**2)
      fftLege = this%fftLege(i:i+step-1)
      
      pmm = 1._dbl
      mj = 0
      
      sumLegendreN = czero
      sumLegendreS = czero

      m = 0
        pmj2 = 0._dbl
        pmj1 = 0._dbl
        pmj  = 1._dbl

        j = m
          s = +1 ; mj = mj + 1

          do concurrent ( i2=1:step, i1=1:6 )
            sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj)
            sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj)
          end do

        do j = m+1, this%maxj
          s = -s ; mj = mj+1
          
          pmj2 = pmj1
          pmj1 = pmj
          pmj  = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2

          do concurrent ( i2=1:step, i1=1:6 )
            sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj) * pmj(i2)
            sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj) * pmj(i2) * s
          end do
        end do

        do concurrent ( i2=1:step, i1=1:6 )
          sumLegendreN(i1,i2,m)%im = 0._dbl
          sumLegendreS(i1,i2,m)%im = 0._dbl
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

            do concurrent ( i2=1:step, i1=1:6 )
              sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj)
              sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj)
            end do

          do j = m+1, this%maxj
            s = -s ; mj = mj+1

            pmj2 = pmj1
            pmj1 = pmj
            pmj  = this%amjrr(mj) * cosx * pmj1 - this%bmjrr(mj) * pmj2

            do concurrent ( i2=1:step, i1=1:6 )
              sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj) * pmj(i2)
              sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj) * pmj(i2) * s
            end do
          end do
        
          do concurrent ( i2=1:step, i1=1:6 )
            sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) * pmm(i2)
            sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) * pmm(i2)
          end do
        end do

      call fftw_execute_dft_c2r(this%fftw_06_c2r, sumLegendreN, grid)

      do concurrent ( i1=0:this%nFourier-1, i2=1:step )
        fft(i2,i1) = grid(1,i2,i1) * grid(4,i2,i1) + &
                   & grid(2,i2,i1) * grid(5,i2,i1) + &
                   & grid(3,i2,i1) * grid(6,i2,i1)
      end do

      call fftw_execute_dft_r2c(this%fftw_01_r2c, fft, fftNC)

      call fftw_execute_dft_c2r(this%fftw_06_c2r, sumLegendreS, grid)

      do concurrent ( i1=0:this%nFourier-1, i2=1:step )
        fft(i2,i1) = grid(1,i2,i1) * grid(4,i2,i1) + &
                   & grid(2,i2,i1) * grid(5,i2,i1) + &
                   & grid(3,i2,i1) * grid(6,i2,i1)
      end do

      call fftw_execute_dft_r2c(this%fftw_01_r2c, fft, fftSC)
      
      pmm = 1._dbl
      mj  = 0

      m = 0
        pmj2 = 0._dbl
        pmj1 = 0._dbl
        pmj  = 1._dbl

        do concurrent ( i2 = 1:step )
          fftNC(i2,m) = fftLege(i2) * fftNC(i2,m) ; fftNC(i2,m)%im = 0._dbl
          fftSC(i2,m) = fftLege(i2) * fftSC(i2,m) ; fftSC(i2,m)%im = 0._dbl
        end do

        j = m
          s = +1 ; mj = mj+1

          cr(mj) = cr(mj) + sum( fftNC(:,m) + fftSC(:,m) )
        
        do j = m+1, this%jmax+1
          s = -s ; mj = mj+1

          pmj2 = pmj1
          pmj1 = pmj
          pmj  = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2

          cr(mj) = cr(mj) + sum( pmj * ( fftNC(:,m) + s * fftSC(:,m) ) )
        end do
        
        do m = 1, this%jmax+1
          fac = sqrt(( 2*m+1._dbl ) / ( 2*m ))
          pmm = -fac * sinx * pmm
          if (maxval(abs(pmm)) < 1.0d-55) exit

          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl

          do concurrent ( i2 = 1:step )
            fftNC(i2,m) = fftLege(i2) * fftNC(i2,m) * pmm(i2)
            fftSC(i2,m) = fftLege(i2) * fftSC(i2,m) * pmm(i2)
          end do

          j = m
            s = +1 ; mj = mj+1

            cr(mj) = cr(mj) + sum( fftNC(:,m) + fftSC(:,m) )
           
          do j = m+1, this%jmax+1
            s = -s ; mj = mj+1

            pmj2 = pmj1
            pmj1 = pmj
            pmj  = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2

            cr(mj) = cr(mj) + sum( pmj * ( fftNC(:,m) + s * fftSC(:,m) ) )
          end do
      end do
    end do

    deallocate( cc, fft, fftNC, fftSC, pmm, pmj, pmj1, pmj2, cosx, sinx, sumLegendreN, sumLegendreS, fftLege, grid )

    fac = 1 / (16 * this%nLegendre**2 * this%nFourier * sqrt(pi) )

    do j = 0, this%jmax
      do m = 0, j
        cjm(j*(j+1)/2+m+1) = cr(m*(this%maxj)-m*(m+1)/2+j+1) * fac
      end do
    end do

    deallocate(cr)

  end function vcvv_fn

  subroutine deallocate_fftw_vcvv_sub(this)
    class(T_lateralGrid), intent(inout) :: this

    call fftw_destroy_plan( this%fftw_06_c2r )
    call fftw_destroy_plan( this%fftw_01_r2c )

  end subroutine deallocate_fftw_vcvv_sub

end submodule vcvv