submodule(Solution) Solution_singles
  implicit none ; contains
  
  module pure complex(kind=dbl) function deviatoric_stress_fn(this, ir, il, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ijm, il
    integer                       :: isp, ist
    
    deviatoric_stress_fn = czero
    
    if ( ijm >= 2 ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      select case (il)
        case (-2)
          if ( this%initsfer ) deviatoric_stress_fn = this%mech(isp+2,ijm)
        case (-1)
          if ( this%inittorr ) deviatoric_stress_fn = this%torr(ist+1,ijm)
        case ( 0)
          if ( this%initsfer ) deviatoric_stress_fn = this%mech(isp+4,ijm)
        case (+1)
          if ( this%inittorr ) deviatoric_stress_fn = this%torr(ist+2,ijm)
        case (+2)
          if ( this%initsfer ) deviatoric_stress_fn = this%mech(isp+5,ijm)
      end select
    end if
    
  end function deviatoric_stress_fn
  
end submodule Solution_singles