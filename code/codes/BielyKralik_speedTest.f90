program BielyKralik_speedTest
  use omp_lib
  use oceanconv
  implicit none

  type(T_oceanConv) :: oceanmodel
  real(kind=dbl)    :: start, end
  
  !Inicializuj vypocet
  call oceanmodel%init_sub()
  
  !Casova slucka
  start = omp_get_wtime()
    call oceanmodel%iter_sub()
  end = omp_get_wtime()

  write(*,*) (end-start) / oceanmodel%n_iter
  
  !Cistenie
  call oceanmodel%deallocate_sub()
  
end program BielyKralik_speedTest
