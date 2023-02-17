submodule (SphericalHarmonics) vcsv
  implicit none

  contains

  subroutine init_vcsv_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer,              allocatable   :: iemb(:), oemb(:)
    real(kind=dbl),       allocatable   :: testField_re(:,:), testField3_re(:,:,:)
    complex(kind=dbl),    allocatable   :: testField(:,:), testField3(:,:,:)

    allocate(testField3(4,step,this%nFourier))
      this%fftw_04_forw = fftw_plan_many_dft( 1, (/this%nFourier/), 4*step, testField3, iemb, 4*step, 1,                &
                                            &                               testField3, oemb, 4*step, 1, +1, fftw_flags )
    deallocate(testField3)

    allocate(testField3_re(step,3,this%nFourier), testField3(step,3,this%nFourier/2+1))
      this%fftw_03_back = fftw_plan_many_dft_r2c( 1, (/this%nFourier/), 3*step, testField3_re, iemb, 3*step, 1,            &
                                                &                               testField3   , oemb, 3*step, 1, fftw_flags )
    deallocate(testField3_re, testField3)

    write(*,*) 'vcsv initialized'

  end subroutine init_vcsv_sub

  function vcsv_fn(this, cajm, cbjml) result(cjml)
    class(T_lateralGrid), intent(in) :: this
    complex(kind=dbl),    intent(in) :: cajm(:), cbjml(:)
    complex(kind=dbl)                :: cjml(this%jmv)
    integer                          :: i, k, j, m, l, mj, lm, jml_int, jm_int, lmj, i1, i2
    real(kind=dbl),      allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), fftLege(:), fft(:,:,:)
    complex(kind=dbl),   allocatable :: cc(:,:), cr(:,:), sumLegendreN(:,:,:), sumLegendreS(:,:,:), fftNC(:,:,:), fftSC(:,:,:)
    
    allocate( cc(4,this%jms2) ); cc = cmplx(0._dbl, 0._dbl, kind=dbl)
    
      mj = 0
        do m = 0, this%maxj
          do j = m, this%maxj
            mj = mj + 1

            do l = abs(j-1), min(this%jmax, j+1)
              lmj= 3*(l*(l+1)/2+m)+j-l

              if (m == 0) then
                cc(2,mj) = cc(2,mj) + conjg( cbjml(lmj+3) ) * cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(l+j)
                cc(3,mj) = cc(3,mj) +        cbjml(lmj+3)   * cleb1_fn(j,m,1,+1,l,m+1)
                cc(4,mj) = cc(4,mj) +        cbjml(lmj  )   * cleb1_fn(j,m,1, 0,l,m  )

              else
                             cc(2,mj) = cc(2,mj) + cbjml(lmj-3) * cleb1_fn(j,m,1,-1,l,m-1)
                if (l > m+0) cc(3,mj) = cc(3,mj) + cbjml(lmj+3) * cleb1_fn(j,m,1,+1,l,m+1)
                if (l > m-1) cc(4,mj) = cc(4,mj) + cbjml(lmj  ) * cleb1_fn(j,m,1, 0,l,m  )
              end if
            end do

            if (j <= this%jmax) cc(1,mj) = cajm( j*(j+1)/2+m+1 )
                                cc(2,mj) =           cc(2,mj) - cc(3,mj)
                                cc(3,mj) = cunit * ( cc(2,mj) + cc(3,mj) + cc(3,mj) )
          end do
        end do
    
      cc(:,1:this%maxj+1) = cc(:,1:this%maxj+1) / 2 ; cc(1,:) = 2 * cc(1,:)

    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), fftLege(step), sumLegendreN(4,step,0:this%nFourier-1), &
            & sumLegendreS(4,step,0:this%nFourier-1), fft(step,3,0:this%nFourier-1), fftNC(step,3,0:this%nFourier/2),          &
            & fftSC(step,3,0:this%nFourier/2), cr(3,this%jms1) ) ; cr = cmplx(0._dbl, 0._dbl, kind=dbl)

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
          
              do concurrent (i2=1:step, i1=1:4)
                sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) + cc(i1,mj) * pmj(i2)
                sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) + cc(i1,mj) * pmj(i2) * (-1)**(j+m)
              end do

              pmj2 = pmj1
              pmj1 = pmj
              pmj  = this%amjrr(mj+1) * cosx * pmj1 - this%bmjrr(mj+1) * pmj2
            end do
        
            do concurrent (i2=1:step, i1=1:4)
              sumLegendreN(i1,i2,m) = sumLegendreN(i1,i2,m) * pmm(i2)
              sumLegendreS(i1,i2,m) = sumLegendreS(i1,i2,m) * pmm(i2)
            end do

            pmm = -sqrt( (1-cosx**2) * (2*m+3) / (m+1) / 2 ) * pmm; if (maxval(abs(pmm)) < 1.0d-55) exit
          end do

        call fftw_execute_dft(this%fftw_04_forw, sumLegendreN, sumLegendreN)
          do concurrent ( k = 0:this%nFourier-1, i2 = 1:step )
            fft(i2,1,k) = sumLegendreN(1,i2,k)%re * sumLegendreN(2,i2,k)%re
            fft(i2,2,k) = sumLegendreN(1,i2,k)%re * sumLegendreN(3,i2,k)%re
            fft(i2,3,k) = sumLegendreN(1,i2,k)%re * sumLegendreN(4,i2,k)%re
          end do
        call fftw_execute_dft_r2c(this%fftw_03_back, fft, fftNC)

        call fftw_execute_dft(this%fftw_04_forw, sumLegendreS, sumLegendreS)
          do concurrent ( k = 0:this%nFourier-1, i2 = 1:step )
            fft(i2,1,k) = sumLegendreS(1,i2,k)%re * sumLegendreS(2,i2,k)%re
            fft(i2,2,k) = sumLegendreS(1,i2,k)%re * sumLegendreS(3,i2,k)%re
            fft(i2,3,k) = sumLegendreS(1,i2,k)%re * sumLegendreS(4,i2,k)%re
          end do
        call fftw_execute_dft_r2c(this%fftw_03_back, fft, fftSC)

        pmm = 1._dbl; mj = 1
          do m = 0, this%jmax+1
            pmj2 = 0._dbl
            pmj1 = 0._dbl
            pmj  = 1._dbl

            do concurrent ( i1=1:3, i2=1:step )
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

            pmm = -sqrt( (1-cosx**2) * (2*m+3) / (m+1) / 2 ) * pmm ; if (maxval(abs(pmm)) < 1.0d-55) exit
          end do
        end do

        cr = cr / 4 / this%nLegendre**2 / this%nFourier / sqrt(pi)

    deallocate( cc, fft, fftNC, fftSC, pmm, pmj, pmj1, pmj2, cosx, sumLegendreN, sumLegendreS )

      do j = 0, this%jmax
        do m = 0, j
          do l = abs(j-1), j+1
            lm = l*(l+1)/2+m+1; jml_int = 3*(j*(j+1)/2+m)+l-j

            if ( m == 0 ) then
              cjml(jml_int) = (       cr(1,lm+1) - cunit * cr(2,lm+1)  ) * cleb1_fn(l,m+1,1,-1,j,m) / 2 + &
                            & (       cr(3,lm  )                       ) * cleb1_fn(l,m+0,1, 0,j,m)     + &
                            & ( conjg(cr(1,lm+1) - cunit * cr(2,lm+1)) ) * cleb1_fn(l,m-1,1,+1,j,m) / 2

            else
              cjml(jml_int) = ( cr(1,lm+1) - cunit * cr(2,lm+1) ) * cleb1_fn(l,m+1,1,-1,j,m) / 2 + &
                            & ( cr(3,lm  )                      ) * cleb1_fn(l,m+0,1, 0,j,m)     - &
                            & ( cr(1,lm-1) + cunit * cr(2,lm-1) ) * cleb1_fn(l,m-1,1,+1,j,m) / 2
            end if
          end do
        end do
      end do

    deallocate(cr)

  end function vcsv_fn

  subroutine deallocate_fftw_vcsv_sub(this)
    class(T_lateralGrid), intent(inout) :: this

    call fftw_destroy_plan( this%fftw_04_forw )
    call fftw_destroy_plan( this%fftw_03_back )

  end subroutine deallocate_fftw_vcsv_sub

end submodule vcsv