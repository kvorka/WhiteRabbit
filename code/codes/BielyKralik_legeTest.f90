program BielyKralik_legeTest
  use omp_lib
  use SphericalHarmonics
  implicit none

  type(T_lateralGrid)            :: sph
  integer                        :: jcut, j, m, jm
  real(kind=dbl)                 :: xrand, tstart, tend
  complex(kind=dbl), allocatable :: scal1(:), scal2(:), scal3(:)
  
  !Inicializuj vypocet
  jcut = 997
  
  call sph%init_sub(jcut)
  
  allocate( scal1(jcut*(jcut+1)/2+jcut+1), scal2(jcut*(jcut+1)/2+jcut+1), scal3(jcut*(jcut+1)/2+jcut+1) )
  
  scal1 = czero
  scal2 = czero ; scal2(1) = cs4pi
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

  tstart = omp_get_wtime()
  call sph%vcsum_sub( scal1, scal2, scal3 )
  tend = omp_get_wtime(); write(*,*) "vcsum time: ", tend-tstart

  write(*,*) 'vcsum: ', maxval( abs( (scal1(:)-scal3(:)) / scal1(:) ) )

  !Cistenie
  call sph%deallocate_sub()

end program BielyKralik_legeTest