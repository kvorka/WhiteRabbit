program BielyKralik_ice
  use IceCrustMod
  implicit none

  type(T_iceCrust)               :: icecr
  integer                        :: jmind
  complex(kind=dbl), allocatable :: flux_up(:)
  
  !Inicializuj vypocet
  call icecr%init_sub(notides=.true.)
  
  !Inicializuj tepelny tok
  allocate( flux_up( icecr%jms ) ) ; flux_up = cmplx(0._dbl, 0._dbl, kind=dbl)
    open(unit=1, file='flux_ocean.dat', status='old', action='read')
      do
        read(1,*) jmind, flux_up(jmind) ; if (jmind == icecr%jms) exit
      end do
    close(1)

  !Casova slucka
  do
    !Spocitaj tvar hranic
    call icecr%time_scheme_sub( flux_up )
  end do
  
  !Cistenie
  deallocate( flux_up )
  call icecr%deallocate_sub()

end program BielyKralik_ice
