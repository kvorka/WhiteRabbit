submodule(Gravity) Gravity_init
  implicit none; contains
  
  module pure subroutine init_gravity_sub(this, gmod, g)
    class(T_gravity), intent(inout) :: this
    character(len=*), intent(in)    :: gmod
    real(kind=dbl),   intent(in)    :: g
    
    this%gmod = gmod ; this%g = g
    
  end subroutine init_gravity_sub
  
  module subroutine set_gravity_sub(this, Dcrust, omega, exc, nlay, subor)
    class(T_gravity),           intent(inout) :: this
    integer,          optional, intent(in)    :: nlay
    real(kind=dbl),   optional, intent(in)    :: Dcrust, omega, exc
    character(len=*), optional, intent(in)    :: subor
    integer                                   :: i

    if (present(Dcrust)) this%Dcrust = Dcrust
    if (present(omega) ) this%omega  = omega
    if (present(exc)   ) this%exc    = exc 

    if ( this%gmod == 'mod' ) then
      this%nlay = nlay
    
      allocate( this%rho(nlay), this%radius(nlay) )
        open(unit=1, file=subor, status='old', action='read')
          do i = this%nlay, 1, -1
            read(1,*) this%rho(i), this%radius(i) 
          end do
        close(1)
    
      this%radius = this%radius * 1e3 / this%Dcrust
    end if
    
  end subroutine set_gravity_sub
  
  module pure subroutine deallocate_gravity_sub(this)
    class(T_gravity), intent(inout) :: this
    
    if ( allocated( this%rho    ) ) deallocate( this%rho    )
    if ( allocated( this%radius ) ) deallocate( this%radius )
    
  end subroutine deallocate_gravity_sub
  
end submodule Gravity_init