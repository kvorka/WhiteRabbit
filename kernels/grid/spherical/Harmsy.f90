module Harmsy
  use Math
  implicit none
  
  integer, parameter, public :: nth = 180
  
  public  :: harmsy_sub, harmsy_Pj0_sub, harmsy_Pj1_sub
  private :: dpmm_sub
  
  contains
  
  subroutine harmsy_sub(jmax_in, spectra_in, data_out)
    integer,            intent(in)  :: jmax_in
    complex(kind=dbl),  intent(in)  :: spectra_in(:)
    real(kind=dbl),     intent(out) :: data_out(:,:)
    integer                         :: it, ip, ij, im
    real(kind=dbl)                  :: costheta
    real(kind=dbl),    allocatable  :: polyLege(:)
    complex(kind=dbl), allocatable  :: sumLege(:), expphi(:), expmul(:)
    
    allocate( sumLege(0:jmax_in), expphi(2*nth), expmul(2*nth), polyLege(jmax_in+1) )
    
    do ip = 1, 2*nth
      expphi(ip) = exp( cunit * (ip-1) * pi / nth )
    end do
    
    do it = 1, nth
      costheta = cos(pi/nth * (it-0.5_dbl)) ; expmul = cone
      
      sumLege(:) = czero
        do im = 0, jmax_in
          call dpmm_sub( costheta, jmax_in, im, polyLege )
          
          do ij = im, jmax_in
            sumLege(im) = sumLege(im) + spectra_in(ij*(ij+1)/2+im+1) * polyLege(ij-im+1)
          end do
        end do
      
      data_out(:,it) = c2r_fn( sumLege(0) )
        do concurrent ( im = 1:jmax_in , ip = 1:2*nth )
          expmul(ip)      = expmul(ip) * expphi(ip)
          data_out(ip,it) = data_out(ip,it) + 2 * c2r_fn( expmul(ip) * sumLege(im) )
        end do
    end do
    
    deallocate( sumLege, expphi, expmul, polyLege )
    
  end subroutine harmsy_sub
  
  subroutine harmsy_Pj0_sub(jmax_in, spectra_in, data_out)
    integer,         intent(in)  :: jmax_in
    real(kind=dbl),  intent(in)  :: spectra_in(:)
    real(kind=dbl),  intent(out) :: data_out(:)
    integer                      :: it
    real(kind=dbl), allocatable  :: polyLege(:)
    
    allocate( polyLege(jmax_in+1) )
    
    do it = 1, nth
      call dpmm_sub( cos(pi/nth * (it-0.5_dbl)), jmax_in, 0, polyLege )
      data_out(it) = sum( spectra_in(1:jmax_in+1) * polyLege(1:jmax_in+1) )
    end do
    
    deallocate( polyLege )

  end subroutine harmsy_Pj0_sub
  
  subroutine harmsy_Pj1_sub(jmax_in, spectra_in, data_out)
    integer,         intent(in)  :: jmax_in
    real(kind=dbl),  intent(in)  :: spectra_in(:)
    real(kind=dbl),  intent(out) :: data_out(:)
    integer                      :: it
    real(kind=dbl), allocatable  :: polyLege(:)
    
    allocate( polyLege(jmax_in+1) )
    
    do it = 1, nth
      call dpmm_sub( cos(pi/nth * (it-0.5_dbl)), jmax_in, 1, polyLege )
      data_out(it) = sum( spectra_in(1:jmax_in) * polyLege(1:jmax_in) )
    end do
    
    deallocate( polyLege )
    
  end subroutine harmsy_Pj1_sub
  
  subroutine dpmm_sub(x, n, m, p)
    real(kind=dbl), intent(in)  :: x
    integer,        intent(in)  :: n, m
    real(kind=dbl), intent(out) :: p(:)
    integer                     :: j
    
    p(1) = 1._dbl
      do j = 1, m
        p(1) = p(1) * (1-x*x) * (2*j+1._dbl) / (2*j)
      end do
    p(1) = ((-1)**m) * sqrt( p(1) / (4*pi) )
    
    if (m < n) then
      p(2) = sqrt(2*m+3._dbl) * x * p(1)
      
      do j = m+2, n
        p(j-m+1) = sqrt( (2*j-1._dbl) * (2*j  +1)           /             (j**2-m**2)   ) * x * p(j-m  ) - &
                 & sqrt( (2*j+1._dbl) * (  j-m-1) * (j+m-1) / ( (2*j-3) * (j**2-m**2) ) )     * p(j-m-1)
      end do
    end if
    
  end subroutine dpmm_sub
  
end module Harmsy