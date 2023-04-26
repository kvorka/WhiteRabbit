submodule(Solution) Solution_spherical
  implicit none

  contains

  pure function temp_jm_fn(this, i) result(temp)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl)             :: temp(this%jms)

    if ( allocated(this%temp) ) then
      temp = this%temp(3*(i-1)+1,:)
    else
      temp = cmplx(0._dbl, 0._dbl, kind=dbl)
    end if

  end function temp_jm_fn

  pure function flux_jml_fn(this, i) result(flux)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl)             :: flux(this%jmv)
    integer                       :: ijm, ijml, ind1

    ind1 = 3*(i-1)+2

    if ( allocated(this%temp) ) then
      ijml = 1
        flux(ijml) = this%temp(ind1+1,1)

      do ijm = 2, this%jms
        ijml = ijml+3

        flux(ijml-2) = this%temp(ind1, ijm)
        flux(ijml-1) = czero
        flux(ijml  ) = this%temp(ind1+1, ijm)
      end do
    else
      flux = czero
    end if

  end function flux_jml_fn

  pure function velocity_jml_fn(this, i) result(velocity)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl)             :: velocity(this%jmv)
    integer                       :: ijm, ijml, sfer_ind1, torr_ind

    ijml = 1
      velocity(ijml) = czero

    sfer_ind1 = 6*(i-1)+1
    torr_ind  = 3*(i-1)+1
      
    if ( allocated(this%mech) .and. allocated(this%torr) ) then
      do ijm = 2, this%jms
        ijml = ijml+3

        velocity(ijml-2) = this%mech(sfer_ind1  , ijm)
        velocity(ijml-1) = this%torr(torr_ind   , ijm)
        velocity(ijml  ) = this%mech(sfer_ind1+1, ijm)
      end do

    else if ( (.not. allocated(this%torr)) .and. allocated(this%mech) ) then
      do ijm = 2, this%jms
        ijml = ijml+3

        velocity(ijml-2) = this%mech(sfer_ind1  , ijm)
        velocity(ijml-1) = czero
        velocity(ijml  ) = this%mech(sfer_ind1+1, ijm)
      end do
      
    else if ( allocated(this%torr) .and. (.not. allocated(this%mech)) ) then
      do ijm = 2, this%jms
        velocity(3*(ijm-1)) = this%torr(torr_ind, ijm)
      end do
    end if

  end function velocity_jml_fn

  pure function conv_velocity_jml_fn(this, i) result(velocity)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl)             :: velocity(this%jmv)
    integer                       :: j, m, l
      
    velocity = cmplx(0._dbl, 0._dbl, kind=dbl)
      
    do j = 0, this%jmax
      do m = 1, j
        do l = abs(j-1)-j, +1
          velocity(3*(j*(j+1)/2+m)+l) = this%velocity_fn(i, j, m, l)
        end do
      end do
    end do

  end function conv_velocity_jml_fn

  pure function deviatoric_stress_jml2_fn(this, i) result(dstress)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl)             :: dstress(this%jmt)
    integer                       :: j, m, l

    do j = 0, this%jmax
      do m = 0, j
        do l = abs(j-2)-j, +2
          dstress(5*(j*(j+1)/2+m)+l-1) = this%deviatoric_stress_fn(i, j, m, l)
        end do
      end do
    end do

  end function deviatoric_stress_jml2_fn

end submodule Solution_spherical