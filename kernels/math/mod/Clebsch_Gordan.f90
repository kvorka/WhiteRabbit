module Clebsch_Gordan
  use Math
  implicit none; public; contains
  
  pure real(kind=dbl) function cleb_fn(j1j, m2, j, m)
    integer, intent(in) :: j1j, m2, j, m
    
    select case (m2)
      case (-1)
        select case (j1j)
          case (-1) ; cleb_fn = +sqrt( (j-m-1) * (j-m  ) / ( (2*j-1._real64) * (  j  ) * 2 ) )
          case ( 0) ; cleb_fn = +sqrt( (j+m+1) * (j-m  ) / ( (  j+1._real64) * (  j  ) * 2 ) )
          case (+1) ; cleb_fn = +sqrt( (j+m+2) * (j+m+1) / ( (  j+1._real64) * (2*j+3) * 2 ) )
        end select
        
      case (0)
        select case (j1j)
          case (-1) ; cleb_fn = +sqrt( (j+m  ) * (j-m  ) / ( (2*j-1._real64) * (  j  ) ) )
          case ( 0) ; cleb_fn = +sqrt( (  m  ) * (  m  ) / ( (  j+1._real64) * (  j  ) ) )
          case (+1) ; cleb_fn = -sqrt( (j+m+1) * (j-m+1) / ( (  j+1._real64) * (2*j+3) ) )
        end select
        
      case (+1)
        select case (j1j)
          case (-1) ; cleb_fn = +sqrt( (j+m-1) * (j+m  ) / ( (2*j-1._real64) * (  j  ) * 2 ) )
          case ( 0) ; cleb_fn = -sqrt( (j+m  ) * (j-m+1) / ( (  j+1._real64) * (  j  ) * 2 ) )
          case (+1) ; cleb_fn = +sqrt( (j-m+1) * (j-m+2) / ( (  j+1._real64) * (2*j+3) * 2 ) )
        end select
    end select
    
  end function cleb_fn
  
end module Clebsch_Gordan
