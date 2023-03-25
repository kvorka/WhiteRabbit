program BielyKralik_speedTest
  use omp_lib
  use OceanConvMod
  implicit none

  type(T_oceanConv) :: ocean
  real(kind=dbl)    :: start, end
  
  !Inicializuj vypocet
  call ocean%init_sub()
  
  !Casova slucka
  start = omp_get_wtime()
    call ocean%speed_sub()
  end = omp_get_wtime()

  write(*,*) (end-start) / ocean%n_iter
  
  !Cistenie
  call ocean%deallocate_sub()

end program BielyKralik_speedTest
