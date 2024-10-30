submodule (radial_grid) interpd
  implicit none ; contains
  
  module procedure c
    
    select case (p)
      case (-1)
        c = ( this%r(i) - this%rr(i+1) ) / ( this%rr(i) - this%rr(i+1) )
      
      case (+1)
        c = ( this%r(i) - this%rr(i) ) / ( this%rr(i+1) - this%rr(i) )
      
    end select
    
  end procedure c
  
  module procedure cc
    
    select case (p)
      case (-1)
        cc = ( this%rr(i) - this%r(i) ) / ( this%r(i-1) - this%r(i) )
      
      case (+1)
        cc = ( this%rr(i) - this%r(i-1) ) / ( this%r(i) - this%r(i-1) )
      
    end select
    
  end procedure cc
  
end submodule interpd