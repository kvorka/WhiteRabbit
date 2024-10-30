submodule (gravity) acceleration
  implicit none; contains
  
  module procedure g_fn
    
    select case(this%gmod)
      case('hom')
        g_fn = this%g

      case('lin')
        g_fn = ri
      
      case('new')
        g_fn = (1/ri)**2
        
    end select

    g_fn = g_fn / this%g
    
  end procedure g_fn
  
end submodule acceleration