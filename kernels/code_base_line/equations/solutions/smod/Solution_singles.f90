submodule(Solution) Solution_singles
  implicit none

  contains

  pure complex(kind=dbl) function temp_fn(this, i, j, m)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: i, j, m

    temp_fn = this%temp(3*(i-1)+1, j*(j+1)/2+m+1)

  end function temp_fn

  pure complex(kind=dbl) function flux_fn(this, i, j, m, l)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: i, j, m, l

    select case (l)
      case (-1)
        flux_fn = this%temp(3*(i-1)+2, j*(j+1)/2+m+1)
      case (0)
        flux_fn = cmplx(0._dbl, 0._dbl, kind=dbl)
      case (+1)
        flux_fn = this%temp(3*(i-1)+3, j*(j+1)/2+m+1)
      end select

  end function flux_fn

  pure complex(kind=dbl) function velocity_fn(this, i, j, m, l)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: i, j, m, l

    velocity_fn = cmplx(0._dbl, 0._dbl, kind=dbl)

    if ( j > 0 ) then
      select case (l)
        case (-1)
          if ( allocated(this%mech) ) velocity_fn = this%mech(6*(i-1)+1, j*(j+1)/2+m+1)
        case ( 0)
          if ( allocated(this%torr) ) velocity_fn = this%torr(3*(i-1)+1, j*(j+1)/2+m+1)
        case (+1)
          if ( allocated(this%mech) ) velocity_fn = this%mech(6*(i-1)+2, j*(j+1)/2+m+1)
      end select
    end if

  end function velocity_fn

  pure complex(kind=dbl) function deviatoric_stress_fn(this, i, j, m, l)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: i, j, m, l

    deviatoric_stress_fn = cmplx(0._dbl, 0._dbl, kind=dbl)

    if ( j > 0 ) then
      select case (l)
        case (-2)
          if ( allocated(this%mech) ) deviatoric_stress_fn = this%mech(6*(i-1)+3, j*(j+1)/2+m+1)
        case (-1)
          if ( allocated(this%torr) ) deviatoric_stress_fn = this%torr(3*(i-1)+2, j*(j+1)/2+m+1)
        case ( 0)
          if ( allocated(this%mech) ) deviatoric_stress_fn = this%mech(6*(i-1)+5, j*(j+1)/2+m+1)
        case (+1)
          if ( allocated(this%torr) ) deviatoric_stress_fn = this%torr(3*(i-1)+3, j*(j+1)/2+m+1)
        case (+2)
          if ( allocated(this%mech) ) deviatoric_stress_fn = this%mech(6*(i-1)+6, j*(j+1)/2+m+1)
        end select
    end if

  end function deviatoric_stress_fn

end submodule Solution_singles