program BielyKralik_ice
  use icecrust
  implicit none

  type(T_iceCrust)               :: icecr
  integer                        :: indx
  real(kind=dbl)                 :: rescale
  complex(kind=dbl)              :: valf
  complex(kind=dbl), allocatable :: flux_up(:)
  
  !! Inicializacia modelu ladovej kory
  call icecr%init_sub()
  
  !! Inicializacia pociatocneho odhadu pre variaciu tepelneho toku
  allocate( flux_up(icecr%jms) )
  
  open( unit=1, file='flux-averaged.spec', status='old', action='read' )
    do
      read(1,*) indx, valf
      
      if ( ( indx == 1 ) .and. ( valf%re > s4pi ) ) then
        rescale = valf%re / s4pi
      else
        rescale = 1._dbl
      end if
      
      if ( indx <= icecr%jms ) then
        flux_up(indx) = valf / rescale
      else
        exit
      end if
    end do
  close(1)
    
  !!Vypocet
  call icecr%solve_sub(flux_up)
  
  !! Nedokonale cistenie po vypocte
  call icecr%deallocate_sub()
  
end program BielyKralik_ice
