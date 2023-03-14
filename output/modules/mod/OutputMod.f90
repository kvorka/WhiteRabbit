module OutputMod
  use Math
  use FFT
  implicit none

  integer,        parameter          :: nth      = 180     !rozlisenie siete
  real(kind=dbl), parameter          :: dtheta   = pi/nth  !krok na sieti
  real(kind=dbl), parameter, private :: degToRad = pi/180  !prevod stupnov na radiany

  public :: toGrid_sub
  public :: out_data_sub
  public :: out_zondata_sub
  public :: out_spectra1_sub
  public :: out_spectra_sub
  public :: zonalVariation_fn
  public :: tangentCylinder_fn

  private :: harmsy_sub
  private :: dpmm_sub

  contains

  subroutine toGrid_sub(jmax, spectra_in, data_out)
    integer,                           intent(in)  :: jmax
    complex(kind=dbl), dimension(:),   intent(in)  :: spectra_in
    real(kind=dbl),    dimension(:,:), intent(out) :: data_out

    call harmsy_sub(jmax, spectra_in, data_out)

  end subroutine toGrid_sub

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

  subroutine out_data_sub(opt, data_in)
    character(len=*), intent(in) :: opt
    real(kind=dbl),   intent(in) :: data_in(:,:)
    integer                      :: i, j
      
    open(unit=1, file=opt, status='new', action='write')
      
    i = nth
      do j = 2, 2*nth
        write(1,*) j-1, +90, data_in(j,i)
      end do
    
    j = 1
      write(1,*) 360, +90, data_in(j,i)
    
    do i = 2, nth
      do j = 2, 2*nth
        write(1,*) j-1, i-91, (data_in(j,i)+data_in(j,i-1))/2
      end do
    
      j = 1
        write(1,*) 360, i-91, (data_in(j,i)+data_in(j,i-1))/2
    end do
    
    i = 1
      do j = 2, 2*nth
        write(1,*) j-1, -90, data_in(j,i)
      end do
    
    j = 1
      write(1,*) 360, -90, data_in(j,i)
      
    close(1)
  
  end subroutine out_data_sub

  subroutine out_zondata_sub(opt, data_in)
    character(len=*), intent(in) :: opt
    real(kind=dbl),   intent(in) :: data_in(:,:)
    integer                      :: i, j
      
    open(unit=1, file=opt, status='new', action='write')
      
    i = 1
      write(1,*) -90, data_in(1,i)
    
    do i = 2, nth
      write(1,*) i-91, (data_in(1,i)+data_in(1,i-1))/2
    end do

    i = nth
      write(1,*) +90, data_in(1,i)
      
    close(1)

  end subroutine out_zondata_sub

  subroutine out_spectra1_sub(opt, data_in)
    character(len=*),  intent(in) :: opt
    complex(kind=dbl), intent(in) :: data_in(:)
    integer                       :: jm

    open(unit=1, file=opt, status='new', action='write')
      do jm = 1, size(data_in)
        write(1,*) jm, data_in(jm)
      end do
    close(1)

  end subroutine out_spectra1_sub

  subroutine out_spectra_sub(opt, r, data_in)
    character(len=*),  intent(in) :: opt
    real(kind=dbl),    intent(in) :: r(:)
    complex(kind=dbl), intent(in) :: data_in(:,:)
    integer                       :: i

    open(unit=1, file=opt, status='new', action='write')
      do i = 1, size(r)
        write(1,*) r(i), data_in(:,i)
      end do
    close(1)

  end subroutine out_spectra_sub

  real(kind=dbl) function zonalVariation_fn(data_in)
    real(kind=dbl), intent(in) :: data_in(:,:)

    zonalVariation_fn = ( maxval(data_in(1,:)) - minval(data_in(1,:)) )/2

  end function zonalVariation_fn

  real(kind=dbl) function tangentCylinder_fn(data_in, tangentTHT)
    real(kind=dbl), intent(in) :: data_in(:,:), tangentTHT
    integer                    :: i
    real(kind=dbl)             :: qh, ql, theta

    qh = 0._dbl; ql = 0._dbl

    do i = 1, nth
      theta = (i-0.5_dbl)*degToRad

      if ((theta < tangentTHT) .or. (theta > (pi-tangentTHT))) then
        qh = qh + (data_in(1,i) + data_in(1,nth-i+1))/2 * sin(theta) * dtheta
      else
        ql = ql + (data_in(1,i) + data_in(1,nth-i+1))/2 * sin(theta) * dtheta
      end if
    end do

    tangentCylinder_fn = (cos(tangentTHT)*qh - (1-cos(tangentTHT))*ql)/(cos(tangentTHT)*qh + (1-cos(tangentTHT))*ql)

  end function tangentCylinder_fn

end module OutputMod