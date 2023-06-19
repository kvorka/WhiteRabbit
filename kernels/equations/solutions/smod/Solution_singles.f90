submodule(Solution) Solution_singles
  implicit none
  
  contains
  
  pure complex(kind=dbl) function temp_fn(this, ir, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ijm
    integer                       :: is
    
    temp_fn = czero
    
    if ( this%inittemp ) then
      is = 3*(ir-1)+1 ; temp_fn = this%temp(is,ijm)
    end if
    
  end function temp_fn
  
  pure complex(kind=dbl) function flux_fn(this, ir, il, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, il, ijm
    integer                       :: is
    
    flux_fn = czero
    
    is = 3*(ir-1)+1
    
    select case (il)
      case (-1)
        if ( this%inittemp ) flux_fn = this%temp(is+1,ijm)
      case ( 0)
        flux_fn = czero
      case (+1)
        if ( this%inittemp ) flux_fn = this%temp(is+2,ijm)
    end select
    
  end function flux_fn
  
  pure complex(kind=dbl) function velocity_fn(this, ir, il, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ijm, il
    integer                       :: isp, ist
    
    velocity_fn = czero
    
    if ( ijm >= 2 ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      select case (il)
        case (-1)
          if ( this%initsfer ) velocity_fn = this%mech(isp,ijm)
        case ( 0)
          if ( this%inittorr ) velocity_fn = this%torr(ist,ijm)
        case (+1)
          if ( this%initsfer ) velocity_fn = this%mech(isp+1,ijm)
      end select
    end if
    
  end function velocity_fn
  
  pure complex(kind=dbl) function deviatoric_stress_fn(this, ir, il, ijm)
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