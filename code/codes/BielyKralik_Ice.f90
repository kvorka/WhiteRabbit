program BielyKralik_ice
  use IceCrustMod
  implicit none

  type(T_iceCrust)               :: icecr
  integer                        :: j, m
  logical                        :: output
  real(kind=dbl)                 :: rvalue, ivalue, rerror, ierror
  complex(kind=dbl), allocatable :: flux_up(:), shape_up(:)
  
  call icecr%init_sub()
  
  allocate( flux_up(icecr%jms) )
  call icecr%solve_sub(flux_up)
  call icecr%deallocate_sub()

end program BielyKralik_ice
