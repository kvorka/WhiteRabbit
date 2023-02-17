submodule (SphericalHarmonics) vcvv
  implicit none

  contains

  subroutine init_vcvv_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer,              allocatable   :: iemb(:), oemb(:)
    real(kind=dbl),       allocatable   :: testField_re(:,:), testField3_re(:,:,:)
    complex(kind=dbl),    allocatable   :: testField(:,:), testField3(:,:,:)

    allocate(testField3(6,step,this%nFourier))
      this%fftw_06_forw = fftw_plan_many_dft( 1, (/this%nFourier/), 6*step, testField3, iemb, 6*step, 1,                &
                                            &                               testField3, oemb, 6*step, 1, +1, fftw_flags )
    deallocate(testField3)

    allocate(testField_re(step,this%nFourier), testField(step,this%nFourier/2+1))
      this%fftw_01_back = fftw_plan_many_dft_r2c( 1, (/this%nFourier/), step, testField_re, iemb, step, 1,            &
                                                &                             testField   , oemb, step, 1, fftw_flags )
    deallocate(testField_re, testField)

    write(*,*) 'vcvv initialized'

  end subroutine init_vcvv_sub

  function vcvv_fn(this, cajml, cbjml) result(cjm)
    class(T_lateralGrid), intent(in) :: this
    complex(kind=dbl),    intent(in) :: cajml(:), cbjml(:)
    complex(kind=dbl)                :: cjm(this%jmax*(this%jmax+1)/2+this%jmax+1)
    integer                          :: i, k, j, m, l, mj, i1, i2
    real(kind=dbl)                   :: cleb1, cleb2, cleb3
    real(kind=dbl),     allocatable  :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), fftLege(:), fft(:,:)
    complex(kind=dbl),  allocatable  :: cc(:,:), fftNC(:,:), fftSC(:,:), sumLegendreN(:,:,:), sumLegendreS(:,:,:)

    cjm = cmplx(0._dbl, 0._dbl, kind=dbl)
    
    allocate( cc(6,this%jms2) ); cc = cmplx(0._dbl, 0._dbl, kind=dbl)
    
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
    
    cc(:,1:this%maxj+1) = cc(:,1:this%maxj+1) / 2
    
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), fftLege(step), sumLegendreN(6,step,0:this%nFourier-1), &
            & sumLegendreS(6,step,0:this%nFourier-1), fft(step,0:this%nFourier-1), fftNC(step,0:this%nFourier/2),              &
            & fftSC(step,0:this%nFourier/2) )

    do i = 1, this%nLegendre, step
      cosx    = this%roots(i:i+step-1)
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

            do concurrent ( i2 = 1:step, i1 = 1:6 )
              sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj) * pmj(i2)
              sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj) * pmj(i2) * (-1)**(j+m)
            end do

            pmj2 = pmj1
            pmj1 = pmj
            pmj  = this%amjrr(mj+1) * cosx * pmj1 - this%bmjrr(mj+1) * pmj2
          end do
        
          do concurrent ( i2 = 1:step, i1 = 1:6 )
            sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) * pmm(i2)
            sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) * pmm(i2)
          end do

          pmm = -sqrt( (1-cosx**2) * (2*m+3) / (m+1) / 2 ) * pmm ; if (maxval(abs(pmm)) < 1.0d-55) exit
        end do

      call fftw_execute_dft(this%fftw_06_forw, sumLegendreN, sumLegendreN)
        do concurrent ( i1 = 0:this%nFourier-1, i2 = 1:step )
          fft(i2,i1) = sumLegendreN(1,i2,i1)%re * sumLegendreN(4,i2,i1)%re + &
                     & sumLegendreN(2,i2,i1)%re * sumLegendreN(5,i2,i1)%re + &
                     & sumLegendreN(3,i2,i1)%re * sumLegendreN(6,i2,i1)%re
        end do
      call fftw_execute_dft_r2c(this%fftw_01_back, fft, fftNC)

      call fftw_execute_dft(this%fftw_06_forw, sumLegendreS, sumLegendreS)
        do concurrent ( i1 = 0:this%nFourier-1, i2 = 1:step )
          fft(i2,i1) = sumLegendreS(1,i2,i1)%re * sumLegendreS(4,i2,i1)%re + &
                     & sumLegendreS(2,i2,i1)%re * sumLegendreS(5,i2,i1)%re + &
                     & sumLegendreS(3,i2,i1)%re * sumLegendreS(6,i2,i1)%re
        end do
      call fftw_execute_dft_r2c(this%fftw_01_back, fft, fftSC)
      
      pmm = 1._dbl; mj = 1
        do m = 0, this%jmax
          pmj2 = 0._dbl
          pmj1 = 0._dbl
          pmj  = 1._dbl

          do concurrent ( i2 = 1:step )
            fftNC(i2,m) = fftLege(i2) * fftNC(i2,m) * pmm(i2)
            fftSC(i2,m) = fftLege(i2) * fftSC(i2,m) * pmm(i2)
          end do

          do j = m, this%jmax
            cjm(j*(j+1)/2+m+1) = cjm(j*(j+1)/2+m+1) + sum( pmj(:) * ( fftNC(:,m) + (-1)**(j+m) * fftSC(:,m) ) )

          mj = mj+1
            pmj2 = pmj1
            pmj1 = pmj
            pmj  = this%amjrr(mj+m) * cosx * pmj1 - this%bmjrr(mj+m) * pmj2
        end do

        mj = mj+1
        pmm = -sqrt( (1-cosx**2) * (2*m+3) / (m+1) / 2 ) * pmm
      end do
    end do

    cjm = cjm / 4 / this%nLegendre**2 / this%nFourier / sqrt(pi)

    deallocate( cc, fft, fftNC, fftSC, pmm, pmj, pmj1, pmj2, cosx, sumLegendreN, sumLegendreS, fftLege )

  end function vcvv_fn

  subroutine deallocate_fftw_vcvv_sub(this)
    class(T_lateralGrid), intent(inout) :: this

    call fftw_destroy_plan( this%fftw_06_forw )
    call fftw_destroy_plan( this%fftw_01_back )

  end subroutine deallocate_fftw_vcvv_sub

end submodule vcvv