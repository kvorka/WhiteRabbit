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
  
  call sph%reindexing%allocate_scalars_sub(1, cin1)
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
  
  call sph%reindexing%allocate_scalars_sub(1, cout)
  call sph%allocate_grid_sub( grid )
    
    call sph%space_to_grid_sub( cin1, grid )
    call sph%grid_to_space_sub( grid, cout )
  
  write(*,*) maxval( abs( cin1 - cout ) )
  
  call sph%deallocate_sub()
  
end program BielyKralik_legeTest