program BielyKralik_oceanConv
    use oceanconv
    implicit none
  
    type(T_oceanConv) :: ocean
    
    !Inicializuj vypocet
    call ocean%init_sub()
    
    !Casova slucka
    do
      call ocean%iter_sub()
    end do
    
    !Cistenie
    call ocean%deallocate_sub()
  
end program BielyKralik_oceanConv
  