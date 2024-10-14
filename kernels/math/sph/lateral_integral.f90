module lateral_integral
  use math
  implicit none; public; contains
  
  pure real(kind=dbl) function theta_average_fn( func )
    interface
      pure real(kind=dbl) function func( theta )
        import                     :: dbl
        real(kind=dbl), intent(in) :: theta
      end function func
    end interface
    
    integer,        parameter :: nt = 1e6
    real(kind=dbl), parameter :: dtheta = pi / nt
    integer        :: it
    real(kind=dbl) :: tsum, tangle
    
    tsum = zero
    
    do it = 0, nt-1
      tangle = (it+half) * dtheta
      tsum  = tsum + sin(tangle) * func(tangle) * dtheta
    end do
    
    theta_average_fn = tsum / 2
    
  end function theta_average_fn
  
end module