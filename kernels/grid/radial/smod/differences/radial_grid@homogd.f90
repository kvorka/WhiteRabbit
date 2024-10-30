submodule (radial_grid) homogd
  implicit none; contains
  
  module procedure hd
    
    select case (p)
      case (-1)
        hd = -1 / ( this%rr(i+1) - this%rr(i) )
        
      case (+1)
        hd = +1 / ( this%rr(i+1) - this%rr(i) )
    end select
    
  end procedure hd
  
  module procedure hdd
    
    select case (p)
      case (-1)
        hdd = -1 / ( this%r(i) - this%r(i-1) )
        
      case (+1)
        hdd = +1 / ( this%r(i) - this%r(i-1) )
    end select
    
  end procedure hdd
  
  module procedure hdrr
    
    select case (p)
      case (-1)
        hdrr = -1 / ( this%rr(i+1) - this%rr(i-1) )
      
      case (0)
        hdrr = zero
      
      case (+1)
        hdrr = +1 / ( this%rr(i+1) - this%rr(i-1) )
        
    end select
    
  end procedure hdrr
  
end submodule homogd