program BielyKralik_legeTest
  use omp_lib
  use lateral_grid
  use sph_norms
  implicit none

  type(T_lateralGrid)            :: sph
  integer                        :: jcut, j, m, l, jm, jml
  real(kind=dbl)                 :: xrand, xrand1, tstart, tend
  complex(kind=dbl), allocatable :: vec1(:), vec2(:), vec3(:), vec4(:), scal1(:), vec1x2(:,:)
  
  !Inicializuj vypocet
  jcut = 213
  
  call sph%init_sub(jcut)
  
  allocate( scal1(jcut*(jcut+1)/2+jcut+1),     &
          & vec1x2(3,jcut*(jcut+1)/2+jcut+1), &
          & vec1(3*(jcut*(jcut+1)/2+jcut)+1),  &
          & vec2(3*(jcut*(jcut+1)/2+jcut)+1),  &
          & vec3(3*(jcut*(jcut+1)/2+jcut)+1),  &
          & vec4(3*(jcut*(jcut+1)/2+jcut)+1)   )
  
  vec1(1) = cmplx(one, zero, kind=dbl)
  vec2(1) = cmplx(-2._dbl, zero, kind=dbl)
  
  do j = 1, jcut/2
    m = 0
      do l = abs(j-1), j+1
        jml = 3*(j*(j+1)/2+m)+l-j
        
        if ( l /= j ) then
          call random_number( xrand ) ; vec1(jml) = cmplx( xrand, zero, kind=dbl)
          call random_number( xrand ) ; vec2(jml) = cmplx( xrand, zero, kind=dbl)
        else
          call random_number( xrand ) ; vec1(jml) = cmplx( zero, xrand, kind=dbl)
          call random_number( xrand ) ; vec2(jml) = cmplx( zero, xrand, kind=dbl)
        end if
      end do
    
    do m = 1, j
      do l = abs(j-1), j+1
        jml = 3*(j*(j+1)/2+m)+l-j
        
        call random_number( xrand ) ; call random_number( xrand1 ) ; vec1(jml) = cmplx( xrand, xrand1, kind=dbl)
        call random_number( xrand ) ; call random_number( xrand1 ) ; vec2(jml) = cmplx( xrand, xrand1, kind=dbl)
      end do
    end do
  end do
  
  call sph%vcvxv_sub(vec1, vec2, vec3)
  call sph%vcvv_sub(vec3, vec2, scal1)
  
  open(unit=1, file='cross_prod.dat', status='old', action='write')
    do jm = 1, jcut*(jcut+1)/2+jcut+1
      write(1,*) jm, scal1(jm), sum( abs(vec1x2(1:3,jm)) ) / 3
    end do
  close(1)
  
  write(*,*) scalnorm2_fn(jcut, scal1)
  
  !Cistenie
  call sph%deallocate_sub()

end program BielyKralik_legeTest