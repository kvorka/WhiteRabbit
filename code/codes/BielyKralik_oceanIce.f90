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
      !Nastav deformaciu oceanskej hranice a iteruj oceanske prudenie
      call ocean%set_boundary_deformation(u_up = icecr%sol%u_dn * icecr%D_ud, t_up = icecr%sol%t_dn * icecr%D_ud)
      call ocean%iter_sub()
      
      !Spocitaj tvar hranic
      call icecr%time_scheme_sub( ocean%flux_up )
    end do
    
    !Cistenie
    call icecr%deallocate_sub()
    call ocean%deallocate_sub()
  
end program BielyKralik_oceanIce
  