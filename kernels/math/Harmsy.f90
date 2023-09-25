module Harmsy
  use Math
  implicit none
  
  integer, parameter, public :: nth = 180
  
  contains
  
  subroutine harmsy_sub(jmax_in, spectra_in, data_out)
    integer,            intent(in)  :: jmax_in
    complex(kind=dbl),  intent(in)  :: spectra_in(:)
    real(kind=dbl),     intent(out) :: data_out(:,:)
    integer                         :: it, ip, j, m
    real(kind=dbl)                  :: costheta
    real(kind=dbl),    allocatable  :: polyLege(:)
    complex(kind=dbl), allocatable  :: sumLege(:), expphi(:), expmul(:)
    
    data_out = 0._dbl
    
    allocate( sumLege(0:jmax_in), expphi(2*nth), expmul(2*nth), polyLege(jmax_in+1) )
    
    do it = 1, nth
      costheta = cos(pi/nth * (it-0.5_dbl)) ; sumLege(:) = czero
      
      do m = 0, jmax_in
        call dpmm_sub( costheta, jmax_in, m, polyLege )
        
        do j = m, jmax_in
          sumLege(m) = sumLege(m) + spectra_in(j*(j+1)/2+m+1) * polyLege(j-m+1)
        end do
      end do
      
      do ip = 1, 2*nth
        expphi(ip) = exp( cunit * (ip-1) * pi / nth )
        expmul(ip) = cone
      end do
      
      m = 0
        data_out(:,it) = data_out(:,it) + c2r_fn( sumLege(m) )
      
      do m = 1, jmax_in
        expmul(:)      = expmul(:) * expphi(:)
        data_out(:,it) = data_out(:,it) + 2 * c2r_fn( expmul(:) * sumLege(m) )
      end do
    end do
      
    deallocate( sumLege, expphi, expmul, polyLege )
      
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