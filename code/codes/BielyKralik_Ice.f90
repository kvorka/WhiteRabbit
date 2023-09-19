program BielyKralik_ice
  use IceCrustMod
  implicit none

  type(T_iceCrust)               :: icecr
  integer                        :: j, m, error
  real(kind=dbl)                 :: rvalue, ivalue, rerror, ierror
  complex(kind=dbl), allocatable :: flux_up(:), shape_up(:)
  
  !! Inicializacia modelu ladovej kory
  call icecr%init_sub()
  
  !! Inicializacia meranych hodnot tvarovej deformacie
  allocate( shape_up(icecr%jms) )
    open(unit=1, file='/mnt/c/Jakub/PhD/WhiteRabbit/output/corlies.cmplx', status='old', action='read')
      do
        read(1,*,iostat=error) j, m, rvalue, ivalue, rerror, ierror
        
        if (error /= 0) then
          exit
        else
          shape_up(jm(j,m)) = cmplx(rvalue, ivalue, kind=dbl) / icecr%D_ud
        end if
      end do
    close(1)
  
  !! Inicializacia pociatocneho odhadu pre variaciu tepelneho toku
  allocate( flux_up(icecr%jms) ) ; flux_up = czero
    flux_up( 4) = cmplx(+0.17_dbl,  0.00_dbl, kind=dbl)
    flux_up( 6) = cmplx(-0.04_dbl, -0.01_dbl, kind=dbl)
    flux_up(11) = cmplx(+0.39_dbl,  0.00_dbl, kind=dbl)
    flux_up(13) = cmplx(+0.06_dbl, -0.12_dbl, kind=dbl)
    flux_up(15) = cmplx(+0.06_dbl, +0.21_dbl, kind=dbl)
  
  !! Inverzia
  do
    call icecr%solve_sub(flux_up) ; write(*,*) icecr%sol%u_up(4) * icecr%D_ud
    
    if ( max( abs( (icecr%sol%u_up( 4)-shape_up( 4)) / shape_up( 4) ) ,        &
       &      abs( (icecr%sol%u_up( 6)-shape_up( 6)) / shape_up( 6) ) ,        &
       &      abs( (icecr%sol%u_up(11)-shape_up(11)) / shape_up(11) ) ,        &
       &      abs( (icecr%sol%u_up(13)-shape_up(13)) / shape_up(13) ) ,        &  
       &      abs( (icecr%sol%u_up(15)-shape_up(15)) / shape_up(15) ) ) < 1e-3 ) exit
    
    flux_up( 4) = ( shape_up( 4) / icecr%sol%u_up( 4) ) * flux_up( 4)
    flux_up( 6) = ( shape_up( 6) / icecr%sol%u_up( 6) ) * flux_up( 6)
    flux_up(11) = ( shape_up(11) / icecr%sol%u_up(11) ) * flux_up(11)
    flux_up(13) = ( shape_up(13) / icecr%sol%u_up(13) ) * flux_up(13)
    flux_up(15) = ( shape_up(15) / icecr%sol%u_up(15) ) * flux_up(15)
  end do

  !! Vypis inverzie
  open(unit=1, file='inversion_output', status='new', action='write')
    write(1,*) 'bottom temperature:'
    write(1,*) icecr%Td
    write(1,*)
    
    write(1,*) 'bottom flux:'
    write(1,*) '0 0: ', c2r_fn( -icecr%sol%flux_fn(1,1,1) / sqrt(4*pi) ) * icecr%lambdaU * (icecr%Td-icecr%Tu) / icecr%D_ud
    write(1,*) '2 0: ', flux_up( 4)
    write(1,*) '2 2: ', flux_up( 6)
    write(1,*) '4 0: ', flux_up(11)
    write(1,*) '4 2: ', flux_up(13)
    write(1,*) '4 4: ', flux_up(15)
    write(1,*)
    
    write(1,*) 'surface shape:'
    write(1,*) '2 0: ', icecr%sol%u_up( 4) * icecr%D_ud
    write(1,*) '2 2: ', icecr%sol%u_up( 6) * icecr%D_ud
    write(1,*) '4 0: ', icecr%sol%u_up(11) * icecr%D_ud
    write(1,*) '4 2: ', icecr%sol%u_up(13) * icecr%D_ud
    write(1,*) '4 4: ', icecr%sol%u_up(15) * icecr%D_ud
    write(1,*)

    write(1,*) 'surface topography:'
    write(1,*) '2 0: ', icecr%sol%t_up( 4) * icecr%D_ud
    write(1,*) '2 2: ', icecr%sol%t_up( 6) * icecr%D_ud
    write(1,*) '4 0: ', icecr%sol%t_up(11) * icecr%D_ud
    write(1,*) '4 2: ', icecr%sol%t_up(13) * icecr%D_ud
    write(1,*) '4 4: ', icecr%sol%t_up(15) * icecr%D_ud
    write(1,*)

    write(1,*) 'bottom topography:'
    write(1,*) '2 0: ', icecr%sol%t_dn( 4) * icecr%D_ud
    write(1,*) '2 2: ', icecr%sol%t_dn( 6) * icecr%D_ud
    write(1,*) '4 0: ', icecr%sol%t_dn(11) * icecr%D_ud
    write(1,*) '4 2: ', icecr%sol%t_dn(13) * icecr%D_ud
    write(1,*) '4 4: ', icecr%sol%t_dn(15) * icecr%D_ud
    write(1,*)

    write(1,*) 'bottom shape:'
    write(1,*) '2 0: ', icecr%sol%u_dn( 4) * icecr%D_ud
    write(1,*) '2 2: ', icecr%sol%u_dn( 6) * icecr%D_ud
    write(1,*) '4 0: ', icecr%sol%u_dn(11) * icecr%D_ud
    write(1,*) '4 2: ', icecr%sol%u_dn(13) * icecr%D_ud
    write(1,*) '4 4: ', icecr%sol%u_dn(15) * icecr%D_ud
    write(1,*)
  close(1)
  
  !! Nedokonale cistenie po vypocte
  call icecr%deallocate_sub()

end program BielyKralik_ice
