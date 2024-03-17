module Clebsch_Gordan
  use Math
  implicit none; public; contains
  
  pure real(kind=dbl) function cleb1_fn(j1, m1, j2, m2, j, m)
    integer, intent(in) :: j1, m1, j2, m2, j, m
    
    select case (m2)
      case (-1)
        select case (j1-j)
          case (-1) ; cleb1_fn = +sqrt( (j-m-1) * (j-m  ) / ( (2*j-one) * (  j  ) * 2 ) )
          case ( 0) ; cleb1_fn = +sqrt( (j+m+1) * (j-m  ) / ( (  j+one) * (  j  ) * 2 ) )
          case (+1) ; cleb1_fn = +sqrt( (j+m+2) * (j+m+1) / ( (  j+one) * (2*j+3) * 2 ) )
        end select
        
      case (0)
        select case (j1-j)
          case (-1) ; cleb1_fn = +sqrt( (j+m  ) * (j-m  ) / ( (2*j-one) * (  j  ) ) )
          case ( 0) ; cleb1_fn = +sqrt( (  m  ) * (  m  ) / ( (  j+one) * (  j  ) ) )
          case (+1) ; cleb1_fn = -sqrt( (j+m+1) * (j-m+1) / ( (  j+one) * (2*j+3) ) )
        end select
        
      case (+1)
        select case (j1-j)
          case (-1) ; cleb1_fn = +sqrt( (j+m-1) * (j+m  ) / ( (2*j-one) * (  j  ) * 2 ) )
          case ( 0) ; cleb1_fn = -sqrt( (j+m  ) * (j-m+1) / ( (  j+one) * (  j  ) * 2 ) )
          case (+1) ; cleb1_fn = +sqrt( (j-m+1) * (j-m+2) / ( (  j+one) * (2*j+3) * 2 ) )
        end select
    end select
    
  end function cleb1_fn
  
  pure real(kind=dbl) function cleb2_fn(j1, m1, j2, m2, j, m)
    integer, intent(in) :: j1, m1, j2, m2, j, m
    integer             :: i
    real(kind=dbl)      :: c2
    
    c2 = zero; cleb2_fn = c2
    
    i  = sign(1,m2) * m
    
    select case ( abs(m2) )
      
      case (0)
        select case (j-j1)
          case (+2); c2 =      sqrt( ( 3 * (j+m-one) * (j+m  ) * (j-m-1) * (j-m  ) ) / ( (2*j-3  ) * (2*j-2) * (2*j-1) * (  j  ) ) )
          case (+1); c2 = +m * sqrt( ( 3 * (j+m    ) * (j-m  )                     ) / ( (  j-one) * (2*j-1) * (  j+1) * (  j  ) ) )
          case (-1); c2 = -m * sqrt( ( 3 * (j+m+one) *           (j-m+1)           ) / ( (2*j+3  ) * (  j+1) * (  j+2) * (  j  ) ) )
          case (-2); c2 =      sqrt( ( 3 * (j+m+one) * (j+m+2) * (j-m+1) * (j-m+2) ) / ( (  j+1  ) * (2*j+3) * (2*j+4) * (2*j+5) ) )
          case ( 0); c2 = ( 3._dbl * m * m - j * (j+1) ) / sqrt( (2*j-1._dbl) * j * (j+1) * (2*j+3) )
        end select
        
      case (1)
        select case (j-j1)
          case(+2); c2 = sqrt( ( (j+i-2._dbl) * (j+i-1) * (j+i) * (j-i) ) / ( ( 2*j-3 ) * (j-1) * (2*j-1) * j ) )
          case(+1); c2 = -m2 * (j-2*i+1) * sqrt( ( (j+i-1._dbl) * (j+i) ) / ( (2*j-2) * (2*j-1) * j * (j+1) ) )
          case( 0); c2 = (1-2*i) * sqrt( ( 3._dbl * (j+i) * (j-i+1) ) / ( 2 * (2*j-1) * j * (j+1) * (2*j+3) ) )
          case(-1); c2 = m2 * (j+2*i) * sqrt( ( (j-i+1._dbl) * (j-i+2) ) / ( j * (j+1) * (2*j+3) * (2*j+4) ) )
          case(-2); c2 = -sqrt( ( (j+i+1._dbl) * (j-i+1) * (j-i+2) * (j-i+3) ) / ( (j+1) * (2*j+3) *(j+2) * (2*j+5) ) )
        end select
      
      case (2)
        select case (j-j1)
          case(+2); c2 = sqrt( ( (j+i-3._dbl) * (j+i-2) * (j+i-1) * (j+i) ) / ( (2*j-3) * (2*j-2) * (2*j-1) * 2 * j) )
          case(+1); c2 = -m2 / abs(m2) * sqrt( ( (j+i-2._dbl) * (j+i-1) * (j+i) * (j-i+1) ) / ( (2*j-2) * (2*j-1) * j * (j+1) ) )
          case( 0); c2 = sqrt( ( 3._dbl * (j+i-1) * (j+i) * (j-i+1) * (j-i+2) ) / ( 2 * (2*j-1) * j * (j+1) * (2*j+3) ) )
          case(-1); c2 = -m2 / abs(m2) * sqrt( ( (j+i) * (j-i+1._dbl) * (j-i+2) * (j-i+3) ) / ( j * (j+1) * (2*j+3) * (2*j+4) ) )
          case(-2); c2 = sqrt( ( (j-i+1._dbl) * (j-i+2) * (j-i+3) * (j-i+4) ) / ( (2*j+2) * (2*j+3) * (2*j+4) * (2*j+5) ) )
        end select
      
    end select
    
    cleb2_fn = c2
    
  end function cleb2_fn
  
end module Clebsch_Gordan
