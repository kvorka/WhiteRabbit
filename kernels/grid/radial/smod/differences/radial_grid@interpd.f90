submodule (radial_grid) interpd
  implicit none ; contains
  
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
  
end submodule interpd