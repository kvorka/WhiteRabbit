submodule(Gravity) Gravity_init
  implicit none; contains
  
  module pure subroutine init_gravity_sub(this, gmod, g)
    class(T_gravity), intent(inout) :: this
    character(len=*), intent(in)    :: gmod
    real(kind=dbl),   intent(in)    :: g
    
    this%gmod = gmod
    this%g    = g
    
  end subroutine init_gravity_sub
  
  module subroutine set_gravity_sub(this, Dcrust, omega, exc)
    class(T_gravity),           intent(inout) :: this
    real(kind=dbl),   optional, intent(in)    :: Dcrust, omega, exc
    
    if ( present(Dcrust) ) this%Dcrust = Dcrust
    if ( present(omega)  ) this%omega  = omega
    if ( present(exc)    ) this%exc    = exc
        
  end subroutine set_gravity_sub
  
  module pure subroutine deallocate_gravity_sub(this)
    class(T_gravity), intent(inout) :: this
  end subroutine deallocate_gravity_sub
  
end submodule Gravity_init