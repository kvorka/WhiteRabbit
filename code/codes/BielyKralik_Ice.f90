program BielyKralik_ice
  use IceCrustMod
  implicit none

  type(T_iceCrust)               :: icecr
  complex(kind=dbl), allocatable :: flux_up(:)
  
  !! Inicializacia modelu ladovej kory
  call icecr%init_sub()
  
  !! Inicializacia pociatocneho odhadu pre variaciu tepelneho toku
  allocate( flux_up(icecr%jms) )
    flux_up    = czero
    flux_up(4) = czero
  
  !!Vypocet
  call icecr%solve_sub(flux_up)
  write(*,*) icecr%bnd%u_up(4) * icecr%D_ud
  
  !! Nedokonale cistenie po vypocte
  call icecr%deallocate_sub()
  
end program BielyKralik_ice
