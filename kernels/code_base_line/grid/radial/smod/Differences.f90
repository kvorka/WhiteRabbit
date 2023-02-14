submodule(RadialGrid) Differences
  implicit none

  contains

  pure real(kind=dbl) function d(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p
    integer                         :: j

    if (i == 1) then
      select case (p)
        case (-2)
          d = 0._dbl

        case (-1)
          d = (this%r(i)-this%rr(i+1)) + (this%r(i)-this%rr(i+2))

          do j = 0, +2
            if (j /= 0) d = d/(this%rr(i)-this%rr(i+j))
          end do

        case (+1)
          d = (this%r(i)-this%rr(i)) + (this%r(i)-this%rr(i+2))

          do j = 0, +2
            if (j /= +1) d = d/(this%rr(i+1)-this%rr(i+j))
          end do

        case (+2)
          d = (this%r(i)-this%rr(i)) + (this%r(i)-this%rr(i+1))

          do j = 0, +2
            if (j /= +2) d = d/(this%rr(i+2)-this%rr(i+j))
          end do

      end select

    else if (i == this%nd) then
      select case (p)
        case (-2)
          d = (this%r(i)-this%rr(i)) + (this%r(i)-this%rr(i+1))

          do j = -1, +1
            if (j /= -1) d = d/(this%rr(i-1)-this%rr(i+j))
          end do

        case (-1)
          d = (this%r(i)-this%rr(i-1)) + (this%r(i)-this%rr(i+1))

          do j = -1, +1
            if (j /= 0) d = d/(this%rr(i)-this%rr(i+j))
          end do

        case (+1)
          d = (this%r(i)-this%rr(i-1)) + (this%r(i)-this%rr(i))

          do j = -1, +1
            if (j /= +1) d = d/(this%rr(i+1)-this%rr(i+j))
          end do

        case (+2)
          d = 0._dbl

      end select

    else
      select case (p)
        case (-2)
          d = ( this%r(i) - this%rr(i+1) )*( this%r(i) - this%rr(i+2) ) + &
            & ( this%r(i) - this%rr(i+2) )*( this%r(i) - this%rr(i)   ) + &
            & ( this%r(i) - this%rr(i)   )*( this%r(i) - this%rr(i+1) )

          do j = -1, +2
            if (j /= -1) d = d/(this%rr(i-1)-this%rr(i+j))
          end do

        case (-1)
          d = ( this%r(i) - this%rr(i+1) )*( this%r(i) - this%rr(i+2) ) + &
            & ( this%r(i) - this%rr(i+2) )*( this%r(i) - this%rr(i-1) ) + &
            & ( this%r(i) - this%rr(i-1) )*( this%r(i) - this%rr(i+1) )

          do j = -1, +2
            if (j /= 0) d = d/(this%rr(i)-this%rr(i+j))
          end do

        case (+1)
          d = ( this%r(i) - this%rr(i)   )*( this%r(i) - this%rr(i+2) ) + &
            & ( this%r(i) - this%rr(i+2) )*( this%r(i) - this%rr(i-1) ) + &
            & ( this%r(i) - this%rr(i-1) )*( this%r(i) - this%rr(i)   )

          do j = -1, +2
            if (j /= +1) d = d/(this%rr(i+1)-this%rr(i+j))
          end do

        case (+2)
          d = ( this%r(i) - this%rr(i)   )*( this%r(i) - this%rr(i+1) ) + &
            & ( this%r(i) - this%rr(i+1) )*( this%r(i) - this%rr(i-1) ) + &
            & ( this%r(i) - this%rr(i-1) )*( this%r(i) - this%rr(i)   )

          do j = -1, +2
            if (j /= +2) d = d/(this%rr(i+2)-this%rr(i+j))
          end do

      end select

    end if

  end function d

  pure real(kind=dbl) function c(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p

    select case (p)
      case (-1)
        c = (this%r(i)-this%rr(i+1))/(this%rr(i)-this%rr(i+1))

      case (+1)
        c = (this%r(i)-this%rr(i))/(this%rr(i+1)-this%rr(i))

    end select

  end function c

  pure real(kind=dbl) function dd(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p
    integer                         :: j

    if (i == 2) then
      select case (p)
        case (-2)
          dd = 0._dbl

        case (-1)
          dd = (this%rr(i)-this%r(i)) + (this%rr(i)-this%r(i+1))

          do j = -1, +1
            if (j /= -1) dd = dd/(this%r(i-1)-this%r(i+j))
          end do

        case (+1)
          dd = (this%rr(i)-this%r(i-1)) + (this%rr(i)-this%r(i+1))

          do j = -1, +1
            if (j /= 0) dd = dd/(this%r(i)-this%r(i+j))
          end do

        case (+2)
          dd = (this%rr(i)-this%r(i-1)) + (this%rr(i)-this%r(i))

          do j = -1, +1
            if (j /= +1) dd = dd/(this%r(i+1)-this%r(i+j))
          end do

      end select

    else if (i == this%nd) then
      select case (p)
        case (-2)
          dd = (this%rr(i)-this%r(i)) + (this%rr(i)-this%r(i-1))

          do j = -2, 0
            if (j /= -2) dd = dd/(this%r(i-2)-this%r(i+j))
          end do

        case (-1)
          dd = (this%rr(i)-this%r(i)) + (this%rr(i)-this%r(i-2))

          do j = -2, 0
            if (j /= -1) dd = dd/(this%r(i-1)-this%r(i+j))
          end do

        case (+1)
          dd = (this%rr(i)-this%r(i-1)) + (this%rr(i)-this%r(i-2))

          do j = -2, 0
            if (j /= 0) dd = dd/(this%r(i)-this%r(i+j))
          end do

        case (+2)
          dd = 0._dbl

      end select

    else
      select case (p)
        case (-2)
          dd = ( this%rr(i) - this%r(i)   )*( this%rr(i) - this%r(i+1) ) + &
             & ( this%rr(i) - this%r(i+1) )*( this%rr(i) - this%r(i-1) ) + &
             & ( this%rr(i) - this%r(i-1) )*( this%rr(i) - this%r(i)   )

          do j = -2, +1
            if (j /= -2) dd = dd/(this%r(i-2)-this%r(i+j))
          end do

        case (-1)
          dd = ( this%rr(i) - this%r(i)   )*( this%rr(i) - this%r(i+1) ) + &
             & ( this%rr(i) - this%r(i+1) )*( this%rr(i) - this%r(i-2) ) + &
             & ( this%rr(i) - this%r(i-2) )*( this%rr(i) - this%r(i)   )

          do j = -2, +1
            if (j /= -1) dd = dd/(this%r(i-1)-this%r(i+j))
          end do

        case (+1)
          dd = ( this%rr(i) - this%r(i-1) )*( this%rr(i) - this%r(i+1) ) + &
             & ( this%rr(i) - this%r(i+1) )*( this%rr(i) - this%r(i-2) ) + &
             & ( this%rr(i) - this%r(i-2) )*( this%rr(i) - this%r(i-1) )

          do j = -2, +1
            if (j /= 0) dd = dd/(this%r(i)-this%r(i+j))
          end do

        case (+2)
          dd = ( this%rr(i) - this%r(i-2) )*( this%rr(i) - this%r(i-1) ) + &
             & ( this%rr(i) - this%r(i-1) )*( this%rr(i) - this%r(i)   ) + &
             & ( this%rr(i) - this%r(i)   )*( this%rr(i) - this%r(i-2) )

          do j = -2, +1
            if (j /= +1) dd = dd/(this%r(i+1)-this%r(i+j))
          end do

      end select

    end if

  end function dd

  pure real(kind=dbl) function cc(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p

    select case (p)
      case (-1)
        cc = (this%rr(i)-this%r(i))/(this%r(i-1)-this%r(i))

      case (+1)
        cc = (this%rr(i)-this%r(i-1))/(this%r(i)-this%r(i-1))

    end select

  end function cc

end submodule Differences