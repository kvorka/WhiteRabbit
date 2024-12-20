submodule (radial_grid) chebyschd
  implicit none; contains
  
  module procedure d
    
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

  end procedure d
  
  module procedure dd
    
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

  end procedure dd
  
  module procedure drr
    
    select case (p)
      case (-1)
        drr = +( this%rr(i+1) - this%rr(i) ) / ( this%rr(i-1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i-1) )
        
      case (0)
        drr = ( this%rr(i-1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i-1) ) - &
            & ( this%rr(i+1) - this%rr(i) ) / ( this%rr(i-1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i-1) )
        
      case (+1)
        drr = -( this%rr(i-1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i-1) )
        
    end select
    
  end procedure drr
  
end submodule chebyschd