submodule (radial_grid) homogd
  implicit none; contains
  
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
  
  module pure real(kind=dbl) function hdrr(this, i, p)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, p
    
    select case (p)
      case (-1)
        hdrr = -1 / ( this%rr(i+1) - this%rr(i-1) )
      
      case (0)
        hdrr = zero
      
      case (+1)
        hdrr = +1 / ( this%rr(i+1) - this%rr(i-1) )
        
    end select
    
  end function hdrr
  
end submodule homogd