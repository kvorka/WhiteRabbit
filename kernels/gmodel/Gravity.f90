module Gravity
  use RadialGrid
  implicit none
  
  type, public :: T_gravity
    
    character(len=3),            private :: gmod
    integer,                     private :: nlay
    real(kind=dbl),              private :: Dcrust, g, exc, omega
    real(kind=dbl), allocatable, private :: rho(:), radius(:)
    
    contains
    
    procedure :: init_sub       => init_gravity_sub
    procedure :: set_sub        => set_gravity_sub
    procedure :: deallocate_sub => deallocate_gravity_sub

    procedure :: g_fn
    procedure :: V_tide_fn
    procedure :: V_bnd_fn
    procedure :: V_rho_fn
    procedure :: V_rt_fn
  
  end type T_gravity
  
  private :: init_gravity_sub
  private :: deallocate_gravity_sub
  
  contains
  
  subroutine init_gravity_sub(this, gmod, g)
    class(T_gravity), intent(inout) :: this
    character(len=*), intent(in)    :: gmod
    real(kind=dbl),   intent(in)    :: g
    
    this%gmod = gmod
    this%g    = g
    
  end subroutine init_gravity_sub
  
  subroutine set_gravity_sub(this, Dcrust, omega, exc, nlay, subor)
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
  
  pure real(kind=dbl) function g_fn(this, ri)
    class(T_gravity), intent(in) :: this
    real(kind=dbl),   intent(in) :: ri
    integer                      :: i, ii
    real(kind=dbl)               :: mass
    
    select case(this%gmod)
      case('hom')
        g_fn = this%g

      case('lin')
        g_fn = ri
      
      case('new')
        g_fn = (1/ri)**2

      case('mod')
        i = 1
        do
          if (ri <= this%radius(i)) then
            exit
          else if (abs(ri-this%radius(i)) < 1.0d-8) then
            exit
          else
            i = i + 1
          end if
        end do
    
        mass = 4 * pi * this%rho(1) * this%radius(1)**3 / 3
          do ii = 2, i-1
            mass = mass + 4 * pi * this%rho(ii) * (this%radius(ii)**3 - this%radius(ii-1)**3) / 3
          end do
    
        if (i /= 1) mass = mass + 4 * pi * this%rho(i) * (ri**3 - this%radius(i-1)**3) / 3
     
        g_fn = kappa * mass / ri**2 * this%Dcrust

    end select

    g_fn = g_fn / this%g
    
  end function g_fn
  
  pure complex(kind=dbl) function V_bnd_fn(this, j, m, ri, rb, drho, ujm)
    class(T_gravity),  intent(in) :: this
    integer,           intent(in) :: j, m
    real(kind=dbl),    intent(in) :: ri, rb, drho
    complex(kind=dbl), intent(in) :: ujm
    
    V_bnd_fn = czero
    
    if (ri <= rb) then 
      V_bnd_fn = 4 * pi * kappa * ri * this%Dcrust * drho * (ri/rb)**(j-1) * ujm / (2*j+1) / this%g    
    else
      V_bnd_fn = 4 * pi * kappa * ri * this%Dcrust * drho * (rb/ri)**(j+2) * ujm / (2*j+1) / this%g
    end if
    
  end function V_bnd_fn
    
  pure complex(kind=dbl) function V_rho_fn(this, j, m, ri, field, rad_grid)
    class(T_gravity),    intent(in) :: this
    integer,             intent(in) :: j, m
    real(kind=dbl),      intent(in) :: ri
    complex(kind=dbl),   intent(in) :: field(:)
    class(T_radialGrid), intent(in) :: rad_grid
    
    V_rho_fn = 4 * pi * kappa * ri * this%Dcrust * rad_grid%intR_fn(field) / (2*j+1) / this%g
  
  end function V_rho_fn  
  
  pure complex(kind=dbl) function V_tide_fn(this, j, m, ri, phase)
    class(T_gravity), intent(in) :: this
    integer,          intent(in) :: j, m
    real(kind=dbl),   intent(in) :: ri, phase
    
    select case (j)
      case (2)
        select case (m)
          case (0)
            V_tide_fn = (this%omega * ri)**2 * this%Dcrust * this%exc / this%g * r2c_fn( -sqrt(9*pi/5)*cos(phase) )
          case (2)
            V_tide_fn = (this%omega * ri)**2 * this%Dcrust * this%exc / this%g * cmplx( +sqrt(27*pi/10)*cos(phase) ,          &
                                                                                      & -sqrt(24*pi/ 5)*sin(phase) , kind=dbl )
          case default
            V_tide_fn = czero
        end select
      
      case default
        V_tide_fn = czero
    end select
        
  end function V_tide_fn
  
  pure complex(kind=dbl) function V_rt_fn(this, j, m, ri)
    class(T_gravity), intent(in) :: this
    integer,          intent(in) :: j, m
    real(kind=dbl),   intent(in) :: ri
    
    select case (j)
      case (2)
        select case (m)
          case (0)
            V_rt_fn = r2c_fn( -( this%omega * ri )**2 * this%Dcrust / this%g * sqrt(5*pi/9 ) )
          case (2)
            V_rt_fn = r2c_fn( +( this%omega * ri )**2 * this%Dcrust / this%g * sqrt(3*pi/10) )
          case default
            V_rt_fn = czero
        end select
      
      case default
        V_rt_fn = czero
    end select
    
  end function V_rt_fn
  
  subroutine deallocate_gravity_sub(this)
    class(T_gravity), intent(inout) :: this
    
    if ( allocated( this%rho    ) ) deallocate( this%rho    )
    if ( allocated( this%radius ) ) deallocate( this%radius )
    
  end subroutine deallocate_gravity_sub

end module Gravity