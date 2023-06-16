submodule(Solution) Solution_singles
  implicit none
  
  contains
  
  pure complex(kind=dbl) function temp_fn(this, ir, ij, im)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ij, im
    integer                       :: is, ijm
    
    is  = 3*(ir-1)+1
    ijm = ij*(ij+1)/2+im+1
    
    temp_fn = this%temp(is,ijm)
    
  end function temp_fn
  
  pure complex(kind=dbl) function temp2_fn(this, ir, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ijm
    integer                       :: is
    
    is = 3*(ir-1)+1
    
    temp2_fn = this%temp(is,ijm)
    
  end function temp2_fn
  
  pure complex(kind=dbl) function flux_fn(this, ir, ij, im, il)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ij, im, il
    integer                       :: is, ijm
    
    is  = 3*(ir-1)+1
    ijm = ij*(ij+1)/2+im+1
    
    select case (il)
      case (-1)
        flux_fn = this%temp(is+1,ijm)
      case (0)
        flux_fn = czero
      case (+1)
        flux_fn = this%temp(is+2,ijm)
    end select
      
  end function flux_fn
  
  pure complex(kind=dbl) function flux2_fn(this, ir, il, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ijm, il
    integer                       :: is
    
    is  = 3*(ir-1)+1
    
    select case (il)
      case (-1)
        flux2_fn = this%temp(is+1,ijm)
      case (0)
        flux2_fn = czero
      case (+1)
        flux2_fn = this%temp(is+2,ijm)
    end select
      
  end function flux2_fn
  
  pure complex(kind=dbl) function velocity_fn(this, ir, ij, im, il)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ij, im, il
    integer                       :: isp, ist, ijm
    
    velocity_fn = czero
    
    if ( ij > 0 ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      ijm = ij*(ij+1)/2+im+1
      
      select case (il)
        case (-1)
          if ( allocated(this%mech) ) velocity_fn = this%mech(isp,ijm)
        case ( 0)
          if ( allocated(this%torr) ) velocity_fn = this%torr(ist,ijm)
        case (+1)
          if ( allocated(this%mech) ) velocity_fn = this%mech(isp+1,ijm)
      end select
    end if
    
  end function velocity_fn
  
  pure complex(kind=dbl) function velocity2_fn(this, ir, il, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ijm, il
    integer                       :: isp, ist
    
    velocity2_fn = czero
    
    if ( ijm > 1 ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      select case (il)
        case (-1)
          if ( allocated(this%mech) ) velocity2_fn = this%mech(isp,ijm)
        case ( 0)
          if ( allocated(this%torr) ) velocity2_fn = this%torr(ist,ijm)
        case (+1)
          if ( allocated(this%mech) ) velocity2_fn = this%mech(isp+1,ijm)
      end select
    end if
    
  end function velocity2_fn
  
  pure complex(kind=dbl) function deviatoric_stress_fn(this, ir, ij, im, il)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ij, im, il
    integer                       :: isp, ist, ijm
    
    deviatoric_stress_fn = czero
    
    if ( ij > 0 ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      ijm = ij*(ij+1)/2+im+1
      
      select case (il)
        case (-2)
          if ( allocated(this%mech) ) deviatoric_stress_fn = this%mech(isp+2,ijm)
        case (-1)
          if ( allocated(this%torr) ) deviatoric_stress_fn = this%torr(ist+1,ijm)
        case ( 0)
          if ( allocated(this%mech) ) deviatoric_stress_fn = this%mech(isp+4,ijm)
        case (+1)
          if ( allocated(this%torr) ) deviatoric_stress_fn = this%torr(ist+2,ijm)
        case (+2)
          if ( allocated(this%mech) ) deviatoric_stress_fn = this%mech(isp+5,ijm)
      end select
    end if
    
  end function deviatoric_stress_fn
  
  pure complex(kind=dbl) function deviatoric_stress2_fn(this, ir, il, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ijm, il
    integer                       :: isp, ist
    
    deviatoric_stress2_fn = czero
    
    if ( ijm > 1 ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      select case (il)
        case (-2)
          if ( allocated(this%mech) ) deviatoric_stress2_fn = this%mech(isp+2,ijm)
        case (-1)
          if ( allocated(this%torr) ) deviatoric_stress2_fn = this%torr(ist+1,ijm)
        case ( 0)
          if ( allocated(this%mech) ) deviatoric_stress2_fn = this%mech(isp+4,ijm)
        case (+1)
          if ( allocated(this%torr) ) deviatoric_stress2_fn = this%torr(ist+2,ijm)
        case (+2)
          if ( allocated(this%mech) ) deviatoric_stress2_fn = this%mech(isp+5,ijm)
      end select
    end if
    
  end function deviatoric_stress2_fn
  
end submodule Solution_singles