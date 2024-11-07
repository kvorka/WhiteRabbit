submodule (cleb) gordan1
  implicit none; contains
  
  module procedure cleb1_fn
    
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
    
  end procedure cleb1_fn
  
end submodule gordan1