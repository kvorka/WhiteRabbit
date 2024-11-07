submodule (cleb) gordan2
  implicit none; contains
  
  module procedure cleb2_fn
    integer             :: i, s
    
    if ( m2 < 0) then
      s = 1
      i = -m
    else if ( m2 > 0) then
      s = -1
      i = +m
    else
      s = 0
      i = 0
    end if
    
    select case ( abs(m2) )
      
      case (0)
        select case (j-j1)
          case (+2); cleb2_fn =      sqrt((3 * (j+m-one) * (j-m-1) * (j+m  ) * (j-m  )) / ((2*j-3) * (2*j-two) * (2*j-1) * (  j  )))
          case (+1); cleb2_fn = +m * sqrt((3 *      one  *           (j+m  ) * (j-m  )) / ((  j-1) * (  j+one) * (2*j-1) * (  j  )))
          case ( 0); cleb2_fn = (3 * m**2 - j * (j+1) ) / sqrt((2*j-one) * j * (j+1) * (2*j+3))
          case (-1); cleb2_fn = -m * sqrt((3 * (j+m+one) *           (j-m+1)          ) / ((2*j+3) * (  j+one) * (  j+2) * (  j  )))
          case (-2); cleb2_fn =      sqrt((3 * (j+m+one) * (j+m+2) * (j-m+1) * (j-m+2)) / ((2*j+3) * (  j+one) * (2*j+4) * (2*j+5)))
        end select
        
      case (1)
        select case (j-j1)
          case(+2); cleb2_fn = +sqrt( ( (j+i-2) * (j+i-one) * (j+i) * (j-i) ) / ( ( 2*j-3 ) * (j-one) * (2*j-1) * j ) )
          case(+1); cleb2_fn = -m2 * (j-2*i+1) * sqrt( (     (j+i  ) * (j+i-one) ) / (     (2*j-two) * j * (j+1) * (2*j-1) ) )
          case( 0); cleb2_fn =       ( -2*i+1) * sqrt( ( 3 * (j+i  ) * (j-i+one) ) / ( 2 * (2*j-one) * j * (j+1) * (2*j+3) ) )
          case(-1); cleb2_fn = +m2 * (j+2*i  ) * sqrt( (     (j-i+2) * (j-i+one) ) / ( 2 * (  j+two) * j * (j+1) * (2*j+3) ) )
          case(-2); cleb2_fn = -sqrt( ( (j-i+2) * (j+i+one) * (j-i+1) * (j-i+3) ) / ( (j+one) * (2*j+3) *(j+2) * (2*j+5) ) )
        end select
      
      case (2)
        select case (j-j1)
          case(+2); cleb2_fn =     sqrt((    (j+i-3) * (j+i-2) * (j+i-one) * (j+i  )) / (4 * (j-one) * (2*j-3) * (2*j-1) * (j  )))
          case(+1); cleb2_fn = s * sqrt((    (j+i-2) * (j-i+1) * (j+i-one) * (j+i  )) / (2 * (j-one) * (  j+1) * (2*j-1) * (j  )))
          case( 0); cleb2_fn =     sqrt((3 * (j-i+1) * (j-i+2) * (j+i-one) * (j+i  )) / (2 * (j+one) * (2*j-1) * (2*j+3) * (j  )))
          case(-1); cleb2_fn = s * sqrt((    (j-i+3) * (j-i+2) * (j-i+one) * (j+i  )) / (2 * (j+two) * (2*j+3) * (  j+1) * (j  )))
          case(-2); cleb2_fn =     sqrt((    (j-i+3) * (j-i+2) * (j-i+one) * (j-i+4)) / (4 * (j+two) * (2*j+3) * (2*j+5) * (j+1)))
        end select
      
    end select
    
  end procedure cleb2_fn
  
end submodule gordan2
