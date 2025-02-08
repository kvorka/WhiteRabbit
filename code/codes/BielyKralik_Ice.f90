program BielyKralik_ice
  use icecrust
  implicit none

  type(T_iceCrust)               :: icecr
  integer                        :: indx, j, m
  real(kind=dbl)                 :: rescale
  complex(kind=dbl)              :: valf
  complex(kind=dbl), allocatable :: flux_up(:)
  
  !! Inicializacia modelu ladovej kory
  call icecr%init_sub()
  
  !! Inicializacia pociatocneho odhadu pre variaciu tepelneho toku
  allocate( flux_up(icecr%jms) ); flux_up = czero
  
  !open( unit=1, file='flux-averaged.spec', status='old', action='read' )
  !  do
  !    read(1,*) indx, valf
  !    
  !    if ( indx <= icecr%jms ) then
  !      flux_up(indx) = valf
  !    else
  !      exit
  !    end if
  !  end do
  !close(1)
  !
  !rescale = flux_up(1)%re / s4pi
  !
  !do j = 0, icecr%jmax
  !  m = 0
  !    flux_up(j*(j+1)/2+m+1) = flux_up(j*(j+1)/2+m+1) / rescale
  !  
  !  do m = 1, j
  !    flux_up(j*(j+1)/2+m+1) = czero
  !  end do
  !end do
  
  !!Vypocet
  call icecr%solve_sub(flux_up)
  
  !! Nedokonale cistenie po vypocte
  call icecr%deallocate_sub()
  
end program BielyKralik_ice
