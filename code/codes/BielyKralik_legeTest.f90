program BielyKralik_legeTest
  use lateral_grid
  use omp_lib
  implicit none
  
  integer, parameter :: jcut = 497, jms = jcut*(jcut+1)/2+jcut+1, nr = 1
  
  integer                        :: ij, im, ijm, ir
  real(kind=dbl)                 :: rval, start, end
  real(kind=dbl),    allocatable :: grid(:,:,:)
  complex(kind=dbl), allocatable :: cin1(:), cin2(:), cout(:)
  type(T_lateralGrid)            :: sph
  
  call sph%init_sub( jcut )
  
  allocate( cin1(jms) )
    do ij = 0, jcut
      do im = 0, ij
        ijm = ij*(ij+1)/2+im+1
        
        if ( im /= 0 ) then
          call random_number( rval ); cin1(ijm)%re = rval
          call random_number( rval ); cin1(ijm)%im = rval
        else
          call random_number( rval ); cin1(ijm)%re = rval
        end if
      end do
    end do
  
  allocate( cin2(jms) ); cin2 = czero
    cin2(1) = cs4pi
  
  allocate( cout(jms) )
    call sph%vcss_sub( cin1, cin2, cout )
  
  write(*,*) maxval( abs( cin1 - cout ) )
  
  call sph%deallocate_sub()
  
end program BielyKralik_legeTest