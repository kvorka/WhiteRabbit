submodule(Solution) Solution_radial

  contains

  pure function temp_i_fn(this, j, m) result(temp)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: j, m
    integer                       :: i, jm
    complex(kind=dbl)             :: temp(this%nd+1)

    jm = j*(j+1)/2+m+1
      do i = 1, this%nd+1
        temp(i) = this%temp(3*(i-1)+1,jm)
      end do

  end function temp_i_fn

  pure function flux_i_fn(this, j, m, l) result(flux)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: j, m, l
    integer                       :: i
    complex(kind=dbl)             :: flux(this%nd)

    do i = 1, this%nd
      flux(i) = this%flux_fn(i, j, m, l)
    end do

  end function flux_i_fn

  pure function velocity_i_fn(this, j, m, l) result(velocity)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: j, m, l
    integer                       :: i
    complex(kind=dbl)             :: velocity(this%nd+1)

    do i = 1, this%nd+1
      velocity(i) = this%velocity_fn(i, j, m, l)
    end do

  end function velocity_i_fn

  pure function deviatoric_stress_i_fn(this, j, m, l) result(dstress)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: j, m, l
    integer                       :: i
    complex(kind=dbl)             :: dstress(this%nd)

    do i = 1, this%nd
      dstress(i) = this%deviatoric_stress_fn(i, j, m, l)
    end do

  end function deviatoric_stress_i_fn

end submodule Solution_radial
