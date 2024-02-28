program BielyKralik_legeTest
  use SphericalHarmonics
  implicit none

  type(T_lateralGrid) :: sph
  
  !Inicializuj vypocet
  call sph%init_sub(105)
  
  !Cistenie
  call sph%deallocate_sub()

end program BielyKralik_legeTest