module cleb
  use math
  implicit none
  
  interface
    module pure real(kind=dbl) function cleb1_fn(j1, m1, j2, m2, j, m)
      integer, intent(in) :: j1, m1, j2, m2, j, m
    end function cleb1_fn
    
    module pure real(kind=dbl) function cleb2_fn(j1, m1, j2, m2, j, m)
      integer, intent(in) :: j1, m1, j2, m2, j, m
    end function cleb2_fn
  end interface
  
end module cleb