submodule (gravity) acceleration
  implicit none; contains
  
  module pure real(kind=dbl) function g_fn(this, ri)
    class(T_gravity), intent(in) :: this
    real(kind=dbl),   intent(in) :: ri
    
    select case(this%gmod)
      case('hom')
        g_fn = this%g

      case('lin')
        g_fn = ri
      
      case('new')
        g_fn = (1/ri)**2
        
    end select

    g_fn = g_fn / this%g
    
  end function g_fn
  
end submodule acceleration