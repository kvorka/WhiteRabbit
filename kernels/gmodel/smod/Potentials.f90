submodule(Gravity) Potentials
  implicit none; contains
  
  module pure complex(kind=dbl) function V_bnd_fn(this, j, m, ri, rb, drho, ujm)
    class(T_gravity),  intent(in) :: this
    integer,           intent(in) :: j, m
    real(kind=dbl),    intent(in) :: ri, rb, drho
    complex(kind=dbl), intent(in) :: ujm
    
    V_bnd_fn = czero
    
    if (ri <= rb) then 
      V_bnd_fn = 4 * pi * kappa * ri * this%Dcrust * drho * (ri/rb)**(j-1) * ujm / (2*j+1) / this%g    
    else
      V_bnd_fn = 4 * pi * kappa * ri * this%Dcrust * drho * (rb/ri)**(j+2) * ujm / (2*j+1) / this%g
    end if
    
  end function V_bnd_fn
  
  module pure complex(kind=dbl) function V_rho_fn(this, j, m, ri, field, rad_grid)
    class(T_gravity),    intent(in) :: this
    integer,             intent(in) :: j, m
    real(kind=dbl),      intent(in) :: ri
    complex(kind=dbl),   intent(in) :: field(:)
    class(T_radialGrid), intent(in) :: rad_grid
    
    V_rho_fn = 4 * pi * kappa * ri * this%Dcrust * rad_grid%intR_fn(field) / (2*j+1) / this%g
  
  end function V_rho_fn  
  
  module pure complex(kind=dbl) function V_tide_fn(this, j, m, ri, phase)
    class(T_gravity), intent(in) :: this
    integer,          intent(in) :: j, m
    real(kind=dbl),   intent(in) :: ri, phase
    
    if ( (j == 2) .and. (m == 0) ) then
      V_tide_fn = (this%omega * ri)**2 * this%Dcrust * this%exc / this%g * r2c_fn( -sqrt(9*pi/5) * cos(phase) )
      
    else if ( (j == 2) .and. (m == 2) ) then
      V_tide_fn = (this%omega * ri)**2 * this%Dcrust * this%exc / this%g * cmplx( +sqrt(27*pi/10) * cos(phase) ,          &
                                                                                & -sqrt(24*pi/ 5) * sin(phase) , kind=dbl )
      
    else
      V_tide_fn = czero
    end if
    
  end function V_tide_fn
  
  module pure complex(kind=dbl) function V_rt_fn(this, j, m, ri)
    class(T_gravity), intent(in) :: this
    integer,          intent(in) :: j, m
    real(kind=dbl),   intent(in) :: ri
    
    if ( (j == 2) .and. (m == 0) ) then
      V_rt_fn = r2c_fn( -( this%omega * ri )**2 * this%Dcrust / this%g * sqrt(5*pi/9 ) )
      
    else if ( (j == 2) .and. (m == 2) ) then
      V_rt_fn = r2c_fn( +( this%omega * ri )**2 * this%Dcrust / this%g * sqrt(3*pi/10) )
      
    else
      V_rt_fn = czero
    end if
    
  end function V_rt_fn
  
end submodule Potentials