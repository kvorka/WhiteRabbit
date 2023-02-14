module ConvergenceCurveMod
  use Math
  implicit none
  
  public :: convergence_curve_ocean_sub
  public :: convergence_curve_ice_sub

  contains
  
  subroutine convergence_curve_ocean_sub(path_out, path_in)
    character(len=*), intent(in) :: path_out, path_in
    integer                      :: n, j, m, error
    real(kind=dbl)               :: L2_up
    complex(kind=dbl)            :: coeff_up

    open(unit=8, file=path_out, status='new', action='write')
      n = 0
      
      do
        open(unit=7, file=path_in//trim(adjustl(int2str_fn(n)))//'.dat', status='old', action='read', iostat=error)
        if (error /= 0) exit
        
          L2_up = 0._dbl

          do
            read(7,*,iostat=error) j, m, coeff_up
            if (error /= 0) exit
            
            if (m == 0) then
              L2_up = L2_up + abs( coeff_up )**2
            else  
              L2_up = L2_up + 2 * abs( coeff_up )**2
            end if
          end do
      
          write(8,*) n, sqrt(L2_up)
        close(7)
        
        n=n+1
      end do
    close(8)

  end subroutine convergence_curve_ocean_sub

  subroutine convergence_curve_ice_sub(path_out, path_in)
    character(len=*), intent(in) :: path_out, path_in
    integer                      :: n, j, m, error
    real(kind=dbl)               :: L2_dn, L2_up, dimTime
    complex(kind=dbl)            :: coeff_dn, coeff_up

    open(unit=8, file=path_out, status='new', action='write'); n = 0
      do
        open(unit=7, file=path_in//trim(adjustl(int2str_fn(n)))//'.dat', status='old', action='read', iostat=error)
        if (error /= 0) exit
    
          L2_dn = 0._dbl; L2_up = 0._dbl

          do
            read(7,*,iostat=error) j, m, coeff_dn, coeff_up
            if (error /= 0) exit
            
            if (m == 0) then
              L2_dn = L2_dn + abs( coeff_dn )**2
              L2_up = L2_up + abs( coeff_up )**2
            else
              L2_dn = L2_dn + 2 * abs( coeff_dn )**2
              L2_up = L2_up + 2 * abs( coeff_up )**2
            end if
          end do
      
          write(8,*) n, sqrt(L2_dn), sqrt(L2_up)
        close(7)
        
        n=n+1
      end do
    close(8)

  end subroutine convergence_curve_ice_sub

end module ConvergenceCurveMod
