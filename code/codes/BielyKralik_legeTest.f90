program BielyKralik_legeTest
  use SphericalHarmonics
  implicit none

  type(T_lateralGrid)            :: sph
  integer                        :: jcut, j, m, jm
  real(kind=dbl)                 :: xrand
  complex(kind=dbl), allocatable :: scal1(:), scal2(:), scal3(:)
  
  !Inicializuj vypocet
  jcut = 497
  
  call sph%init_sub(jcut)
  
  allocate( scal1(jcut*(jcut+1)/2+jcut+1), scal2(jcut*(jcut+1)/2+jcut+1), scal3(jcut*(jcut+1)/2+jcut+1) )
  
  scal1 = czero
  scal2 = czero ; scal2(1) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
  scal3 = czero
  
  do j = 0, jcut
    do m = 0, j
      jm = j*(j+1)/2+m+1
      
      if ( m /= 0 ) then
        call random_number( xrand ) ; scal1(jm)%re = xrand
        call random_number( xrand ) ; scal1(jm)%im = xrand
      else
        call random_number( xrand ) ; scal1(jm)%re = xrand
        scal1(jm)%im = zero
      end if
    end do
  end do

  call sph%vcsum_sub( scal1, scal2, scal3 )

  write(*,*) 'vcsum: ', maxval(abs(scal1(:)-scal3(:)))
  
  !Cistenie
  call sph%deallocate_sub()

end program BielyKralik_legeTest