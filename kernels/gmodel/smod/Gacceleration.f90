submodule(Gravity) Gacceleration
  implicit none; contains
  
  module pure real(kind=dbl) function g_fn(this, ri)
    class(T_gravity), intent(in) :: this
    real(kind=dbl),   intent(in) :: ri
    integer                      :: i, ii
    real(kind=dbl)               :: mass
    
    select case(this%gmod)
      case('hom')
        g_fn = this%g

      case('lin')
        g_fn = ri
      
      case('new')
        g_fn = (1/ri)**2

      case('mod')
        i = 1
        do
          if (ri <= this%radius(i)) then
            exit
          else if (abs(ri-this%radius(i)) < 1.0d-8) then
            exit
          else
            i = i + 1
          end if
        end do
    
        mass = 4 * pi * this%rho(1) * this%radius(1)**3 / 3
          do ii = 2, i-1
            mass = mass + 4 * pi * this%rho(ii) * (this%radius(ii)**3 - this%radius(ii-1)**3) / 3
          end do
    
        if (i /= 1) mass = mass + 4 * pi * this%rho(i) * (ri**3 - this%radius(i-1)**3) / 3
     
        g_fn = kappa * mass / ri**2 * this%Dcrust

    end select

    g_fn = g_fn / this%g
    
  end function g_fn
  
end submodule Gacceleration