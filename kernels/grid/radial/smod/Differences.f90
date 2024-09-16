submodule (RadialGrid) Differences
  implicit none ; contains
  
  module pure real(kind=dbl) function hd(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p
    
    select case (p)
      case (-1)
        hd = -1 / ( this%rr(i+1) - this%rr(i) )
        
      case (+1)
        hd = +1 / ( this%rr(i+1) - this%rr(i) )
    end select
    
  end function hd
  
  module pure real(kind=dbl) function d(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p
    
    if ( i == 1 ) then
      select case (p)
        case (-2)
          d = zero
        
        case (-1)
          d = ( ( this%r(i) - this%rr(i+1) ) +                                   &
              & ( this%r(i) - this%rr(i+2) ) ) / ( this%rr(i  ) - this%rr(i+1) ) &
              &                                / ( this%rr(i  ) - this%rr(i+2) )
        
        case (+1)
          d = ( ( this%r(i) - this%rr(i)   ) +                                   &
              & ( this%r(i) - this%rr(i+2) ) ) / ( this%rr(i+1) - this%rr(i)   ) &
              &                                / ( this%rr(i+1) - this%rr(i+2) )
        
        case (+2)
          d = ( ( this%r(i) - this%rr(i)   ) +                                   &
              & ( this%r(i) - this%rr(i+1) ) ) / ( this%rr(i+2) - this%rr(i)   ) &
              &                                / ( this%rr(i+2) - this%rr(i+1) )
        
      end select
    
    else if ( i == this%nd ) then
      select case (p)
        case (-2)
          d = ( ( this%r(i) - this%rr(i)   ) +                                   &
              & ( this%r(i) - this%rr(i+1) ) ) / ( this%rr(i-1) - this%rr(i  ) ) &
              &                                / ( this%rr(i-1) - this%rr(i+1) )
          
        case (-1)
          d = ( ( this%r(i) - this%rr(i-1) ) +                                   &
              & ( this%r(i) - this%rr(i+1) ) ) / ( this%rr(i  ) - this%rr(i-1) ) &
              &                                / ( this%rr(i  ) - this%rr(i+1) )
          
        case (+1)
          d = ( ( this%r(i) - this%rr(i-1) ) +                                   &
              & ( this%r(i) - this%rr(i)   ) ) / ( this%rr(i+1) - this%rr(i-1) ) &
              &                                / ( this%rr(i+1) - this%rr(i  ) )
          
        case (+2)
          d = zero
        
      end select
      
    else
      select case (p)
        case (-2)
          d = ( ( this%r(i) - this%rr(i+1) ) * ( this%r(i) - this%rr(i+2) ) +                                   &
              & ( this%r(i) - this%rr(i+2) ) * ( this%r(i) - this%rr(i  ) ) +                                   &
              & ( this%r(i) - this%rr(i  ) ) * ( this%r(i) - this%rr(i+1) ) ) / ( this%rr(i-1) - this%rr(i  ) ) &
              &                                                               / ( this%rr(i-1) - this%rr(i+1) ) &
              &                                                               / ( this%rr(i-1) - this%rr(i+2) )
          
        case (-1)
          d = ( ( this%r(i) - this%rr(i+1) ) * ( this%r(i) - this%rr(i+2) ) +                                   &
              & ( this%r(i) - this%rr(i+2) ) * ( this%r(i) - this%rr(i-1) ) +                                   &
              & ( this%r(i) - this%rr(i-1) ) * ( this%r(i) - this%rr(i+1) ) ) / ( this%rr(i  ) - this%rr(i-1) ) &
              &                                                               / ( this%rr(i  ) - this%rr(i+1) ) &
              &                                                               / ( this%rr(i  ) - this%rr(i+2) )
          
        case (+1)
          d = ( ( this%r(i) - this%rr(i)   ) * ( this%r(i) - this%rr(i+2) ) +                                   &
              & ( this%r(i) - this%rr(i+2) ) * ( this%r(i) - this%rr(i-1) ) +                                   &
              & ( this%r(i) - this%rr(i-1) ) * ( this%r(i) - this%rr(i)   ) ) / ( this%rr(i+1) - this%rr(i-1) ) &
              &                                                               / ( this%rr(i+1) - this%rr(i  ) ) &
              &                                                               / ( this%rr(i+1) - this%rr(i+2) )
          
        case (+2)
          d = ( ( this%r(i) - this%rr(i)   ) * ( this%r(i) - this%rr(i+1) ) +                                   &
              & ( this%r(i) - this%rr(i+1) ) * ( this%r(i) - this%rr(i-1) ) +                                   &
              & ( this%r(i) - this%rr(i-1) ) * ( this%r(i) - this%rr(i)   ) ) / ( this%rr(i+2) - this%rr(i-1) ) &
              &                                                               / ( this%rr(i+2) - this%rr(i  ) ) &
              &                                                               / ( this%rr(i+2) - this%rr(i+1) )
        
      end select
    end if

  end function d
  
  module pure real(kind=dbl) function c(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p
    
    select case (p)
      case (-1)
        c = ( this%r(i) - this%rr(i+1) ) / ( this%rr(i) - this%rr(i+1) )
      
      case (+1)
        c = ( this%r(i) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i) )
      
    end select
    
  end function c
  
  module pure real(kind=dbl) function hdd(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p
    
    select case (p)
      case (-1)
        hdd = -1 / ( this%r(i) - this%r(i-1) )
        
      case (+1)
        hdd = +1 / ( this%r(i) - this%r(i-1) )
    end select
    
  end function hdd
  
  module pure real(kind=dbl) function dd(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p
    
    if ( i == 2 ) then
      select case (p)
        case (-2)
          dd = zero
        
        case (-1)
          dd = ( ( this%rr(i) - this%r(i)   ) +                                 &
               & ( this%rr(i) - this%r(i+1) ) ) / ( this%r(i-1) - this%r(i  ) ) &
               &                                / ( this%r(i-1) - this%r(i+1) )
        
        case (+1)
          dd = ( ( this%rr(i) - this%r(i-1) ) +                                 &
               & ( this%rr(i) - this%r(i+1) ) ) / ( this%r(i  ) - this%r(i-1) ) &
               &                                / ( this%r(i  ) - this%r(i+1) )
        
        case (+2)
          dd = ( ( this%rr(i) - this%r(i-1) ) +                                 &
               & ( this%rr(i) - this%r(i)   ) ) / ( this%r(i+1) - this%r(i-1) ) &
               &                                / ( this%r(i+1) - this%r(i  ) )
        
      end select
    
    else if ( i == this%nd ) then
      select case (p)
        case (-2)
          dd = ( ( this%rr(i) - this%r(i)   ) +                                 &
               & ( this%rr(i) - this%r(i-1) ) ) / ( this%r(i-2) - this%r(i-1) ) &
               &                                / ( this%r(i-2) - this%r(i  ) )
        
        case (-1)
          dd = ( ( this%rr(i) - this%r(i)   ) +                                 &
               & ( this%rr(i) - this%r(i-2) ) ) / ( this%r(i-1) - this%r(i-2) ) &
               &                                / ( this%r(i-1) - this%r(i  ) )
        
        case (+1)
          dd = ( ( this%rr(i) - this%r(i-1) ) +                               &
               & ( this%rr(i) - this%r(i-2) ) ) / ( this%r(i) - this%r(i-2) ) &
               &                                / ( this%r(i) - this%r(i-1) )
        
        case (+2)
          dd = zero
        
      end select
    
    else
      select case (p)
        case (-2)
          dd = ( ( this%rr(i) - this%r(i  ) ) * ( this%rr(i) - this%r(i+1) ) +                                 &
               & ( this%rr(i) - this%r(i+1) ) * ( this%rr(i) - this%r(i-1) ) +                                 &
               & ( this%rr(i) - this%r(i-1) ) * ( this%rr(i) - this%r(i  ) ) ) / ( this%r(i-2) - this%r(i-1) ) &
               &                                                               / ( this%r(i-2) - this%r(i  ) ) &
               &                                                               / ( this%r(i-2) - this%r(i+1) )
        
        case (-1)
          dd = ( ( this%rr(i) - this%r(i  ) ) * ( this%rr(i) - this%r(i+1) ) +                                 &
               & ( this%rr(i) - this%r(i+1) ) * ( this%rr(i) - this%r(i-2) ) +                                 &
               & ( this%rr(i) - this%r(i-2) ) * ( this%rr(i) - this%r(i  ) ) ) / ( this%r(i-1) - this%r(i-2) ) &
               &                                                               / ( this%r(i-1) - this%r(i  ) ) &
               &                                                               / ( this%r(i-1) - this%r(i+1) )
        
        case (+1)
          dd = ( ( this%rr(i) - this%r(i-1) ) * ( this%rr(i) - this%r(i+1) ) +                               &
               & ( this%rr(i) - this%r(i+1) ) * ( this%rr(i) - this%r(i-2) ) +                               &
               & ( this%rr(i) - this%r(i-2) ) * ( this%rr(i) - this%r(i-1) ) ) / ( this%r(i) - this%r(i-2) ) &
               &                                                               / ( this%r(i) - this%r(i-1) ) &
               &                                                               / ( this%r(i) - this%r(i+1) )
        
        case (+2)
          dd = ( ( this%rr(i) - this%r(i-2) ) * ( this%rr(i) - this%r(i-1) ) +                                 &
               & ( this%rr(i) - this%r(i-1) ) * ( this%rr(i) - this%r(i  ) ) +                                 &
               & ( this%rr(i) - this%r(i  ) ) * ( this%rr(i) - this%r(i-2) ) ) / ( this%r(i+1) - this%r(i-2) ) &
               &                                                               / ( this%r(i+1) - this%r(i-1) ) &
               &                                                               / ( this%r(i+1) - this%r(i  ) )
        
      end select
    end if

  end function dd
  
  module pure real(kind=dbl) function cc(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p
    
    select case (p)
      case (-1)
        cc = ( this%rr(i) - this%r(i) ) / ( this%r(i-1) - this%r(i) )
      
      case (+1)
        cc = ( this%rr(i) - this%r(i-1) ) / ( this%r(i) - this%r(i-1) )
      
    end select
    
  end function cc
  
  module pure real(kind=dbl) function drr(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p
    
    select case (p)
      case (-1)
        drr = +( this%rr(i+1) - this%rr(i) ) / ( this%rr(i-1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i-1) )
        
      case (0)
        drr = ( this%rr(i-1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i-1) ) - &
            & ( this%rr(i+1) - this%rr(i) ) / ( this%rr(i-1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i-1) )
        
      case (+1)
        drr = -( this%rr(i-1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i-1) )
        
    end select
    
  end function drr
  
end submodule Differences