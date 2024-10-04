program BielyKralik_tidesOcean
  use oceantides
  implicit none

  type(T_oceanTides) :: ocean
  
  !Inicializuj vypocet
  call ocean%init_sub()

  !Casova slucka
  do
    call ocean%iter_sub()
  end do

  !Cistenie
  call ocean%deallocate_sub()

end program BielyKralik_tidesOcean
