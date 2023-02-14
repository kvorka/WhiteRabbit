program BielyKralik_oceanIce
    use OceanIceMod
    use IceCrustMod
    implicit none
  
    type(T_iceCrust) :: icecr
    type(T_oceanice) :: ocean
    
    !Inicializuj vypocet
    call icecr%init_sub(notides=.true.)
    call ocean%init_sub()
    
    !Casova slucka
    do
      !Iteruj oceanske prudenie
      call ocean%iter_sub( t_bnd = icecr%sol%t_dn * icecr%D_ud, &
                         & u_bnd = icecr%sol%u_dn * icecr%D_ud  )
      
      !Spocitaj tvar hranic
      call icecr%time_scheme_sub( ocean%flux_up )
    end do
    
    !Cistenie
    call icecr%deallocate_sub()
    call ocean%deallocate_sub()
  
end program BielyKralik_oceanIce
  