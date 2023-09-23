module Harmsy
  use FFT
  implicit none
  
  integer, parameter, public :: nth = 180
  
  contains
  
  subroutine harmsy_sub(jmax, coef, data_out)
    integer,                              intent(in)  :: jmax
    complex(kind=dbl), dimension(:),      intent(in)  :: coef
    real(kind=dbl),    dimension(:,:),    intent(out) :: data_out
    integer                                           :: i, j, m, zn
    real(kind=dbl)                                    :: xth, pnm(jmax+1), re(2*nth+1,2), im(2*nth+1,2)
    complex(kind=dbl)                                 :: csum(2)
  
    call init_fft_sub(2*nth)

    do i = 1, nth/2
      re = 0._dbl; im = 0._dbl; xth = cos(2*(2*(nth-i)+1)*atan(1._dbl)/nth)

      do m = 0, jmax
        csum = cmplx(0._dbl, 0._dbl, kind=dbl)
        call dpmm_sub(xth, jmax, m, pnm)
        
        do j = m, jmax
          csum(1) = csum(1) + coef(j*(j+1)/2+m+1) * pnm(j-m+1)
          csum(2) = csum(2) + coef(j*(j+1)/2+m+1) * pnm(j-m+1) * (-1)**(j+m)
        end do

        if (m == 0) then
          re(m+1,:) = +csum%re / 2; im(m+1,:) = 0._dbl
        else
          re(m+1,:) = +csum%re    ; im(m+1,:) = -csum%im
        end if
      end do

      call fft_sub(re(:,1), im(:,1))
      call fft_sub(re(:,2), im(:,2))
        
      data_out(:,    i  ) = 2 * re(1:2*nth,1) * sqrt( real(2*nth, kind=dbl) )
      data_out(:,nth-i+1) = 2 * re(1:2*nth,2) * sqrt( real(2*nth, kind=dbl) )
    end do
  
  end subroutine harmsy_sub
  
  subroutine dpmm_sub(x, n, m, p)
    real(kind=dbl),               intent(in)  :: x
    integer,                      intent(in)  :: n, m
    real(kind=dbl), dimension(:), intent(out) :: p
    integer                                   :: mp, j, i
    real(kind=dbl)                            :: sth, sou, sthm
  
    sth = sqrt(1-x*x)
  
    sou  = 1._dbl
    sthm = 1._dbl
      do mp = 1, m
        sou = sou*(2*mp+1._dbl)/(2*mp)
        sthm = sthm*sth
  
        if (sthm < 1.0d-55) then
          sthm = 0._dbl
          exit
        end if
      end do
  
    p(1) = ((-1)**m)*sqrt(sou/(4*pi))*sthm
  
    if (m < n) p(2) = sqrt(2*m+3._dbl)*x*p(1)
  
    i = 2
    do j = m+2, n
      i = i+1
      p(i)= sqrt((2*j-1._dbl)*(2*j  +1._dbl)                          /(j-m)/(j+m))*x*p(i-1) - &
          & sqrt((2*j+1._dbl)*(  j-m-1._dbl)*(j+m-1._dbl)/(2*j-3._dbl)/(j-m)/(j+m))  *p(i-2)
    end do
  
  end subroutine dpmm_sub
  
end module Harmsy