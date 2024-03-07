submodule (Fourier_transform) fxtal
  implicit none; contains
  
  module pure recursive subroutine fxztal(m, k, l, x, t, ic, itsum, is, j1, it1)
    integer,        intent(in)    :: m, k, l, ic, itsum, is, j1, it1
    real(kind=dbl), intent(in)    :: t(2,0:*)
    real(kind=dbl), intent(inout) :: x(m,2,0:*)
    integer                       :: ip, j, isd, kd, ld, k1, jd, icd, icdd, it1d
    
    kd = k ; ld = l ; k1 = 1 ; j = j1 ; it1d = it1 ; ip = mod(it1d,4)+2
    
    if ( ic /= 1 ) then
      select case (ip)
        case (4)
          call fxzm4a(m, k1, ld, x(1,1,j*ld), t(1,is+j*(ip-1)))
        case (2)
          call fxzm2a(m, k1, ld, x(1,1,j*ld), t(1,is+j*(ip-1)))
        case (3)
          call fxzm3a(m, k1, ld, x(1,1,j*ld), t(1,is+j*(ip-1)))
        case (5)
          call fxzm5a(m, k1, ld, x(1,1,j*ld), t(1,is+j*(ip-1)))
      end select
    end if
    
    k1 = k1 * ip ; ld = ld / ip ; isd = is + kd * (ip-1) ; kd = kd * ip ; j = j * ip ; it1d = it1d / 4 ; icd = ic + 1
    
    if ( (m * ip * ld > 2048) .and. (icd <= itsum) ) then
      do jd = 0, ip-1
         call fxztal(m, kd, ld, x, t, icd, itsum, isd, j+jd, it1d)
      end do
    else
      do icdd = icd, itsum
        ip = mod(it1d,4)+2
        
        select case (ip)
          case (4)
            call fxzm4a(m, k1, ld, x(1,1,j*ld), t(1,isd+j*(ip-1)))
          case (2)
            call fxzm2a(m, k1, ld, x(1,1,j*ld), t(1,isd+j*(ip-1)))
          case (3)
            call fxzm3a(m, k1, ld, x(1,1,j*ld), t(1,isd+j*(ip-1)))
          case (5)
            call fxzm5a(m, k1, ld, x(1,1,j*ld), t(1,isd+j*(ip-1)))
        end select
        
        k1 = k1 * ip ; ld = ld/ip ; isd = isd + kd * (ip-1) ; kd = kd * ip ; j = j * ip ; it1d = it1d / 4
      end do
    end if
    
  end subroutine fxztal
  
end submodule fxtal