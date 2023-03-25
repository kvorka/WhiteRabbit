module RadialGrid
  use Math
  implicit none

  type, public :: T_radialGrid
    integer                     :: nd
    real(kind=dbl)              :: volume
    real(kind=dbl), allocatable :: r(:), rr(:)
    
    contains
    
    procedure, pass :: init_sub       => init_grid_sub
    procedure, pass :: deallocate_sub => deallocate_grid_sub
    
    procedure, pass :: d 
    procedure, pass :: dd
    procedure, pass :: c 
    procedure, pass :: cc
    procedure, pass :: drr
    procedure, pass :: interpolation_fn            
    procedure, pass :: volumetric_integral_real_fn 
    procedure, pass :: volumetric_integral_cmplx_fn
    procedure, pass :: radial_integral_real_fn
    procedure, pass :: radial_integral_cmplx_fn
    
    generic :: intV_fn => volumetric_integral_real_fn, volumetric_integral_cmplx_fn
    generic :: intR_fn => radial_integral_real_fn    , radial_integral_cmplx_fn
    
  end type T_radialGrid

  interface
    module pure real(kind=dbl) function d(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function d

    module pure real(kind=dbl) function c(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function c

    module pure real(kind=dbl) function dd(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function dd

    module pure real(kind=dbl) function cc(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function cc

    module pure real(kind=dbl) function drr(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function drr
  end interface
  
  interface
    module pure real(kind=dbl) function radial_integral_real_fn(this, field)
      class(T_radialGrid), intent(in) :: this
      real(kind=dbl),      intent(in) :: field(:)
    end function radial_integral_real_fn

    module pure real(kind=dbl) function volumetric_integral_real_fn(this, field)
      class(T_radialGrid), intent(in) :: this
      real(kind=dbl),      intent(in) :: field(:)
    end function volumetric_integral_real_fn

    module pure complex(kind=dbl) function radial_integral_cmplx_fn(this, field)
      class(T_radialGrid), intent(in) :: this
      complex(kind=dbl),   intent(in) :: field(:)
    end function radial_integral_cmplx_fn

    module pure complex(kind=dbl) function volumetric_integral_cmplx_fn(this, field)
      class(T_radialGrid), intent(in) :: this
      complex(kind=dbl),   intent(in) :: field(:)
    end function volumetric_integral_cmplx_fn
  end interface

  interface
    module pure function interpolation_fn(this, dimOut, i, rr1, field) result(resField)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, dimOut
      real(kind=dbl),      intent(in) :: rr1(:)
      complex(kind=dbl),   intent(in) :: field(:,:)
      complex(kind=dbl)               :: resField(dimOut)
    end function interpolation_fn
  end interface

  private :: init_grid_sub
  private :: deallocate_grid_sub

  contains

  subroutine init_grid_sub(this, nr, rd, ru, grid_type)
    class(T_radialGrid), intent(inout) :: this
    integer,             intent(in)    :: nr
    real(kind=dbl),      intent(in)    :: rd, ru
    character(len=5),    intent(in)    :: grid_type
    integer                            :: i
    real(kind=dbl)                     :: dr

    this%nd = nr; allocate( this%r(this%nd), this%rr(this%nd+1) )
    
    select case (grid_type)
      case('chebv')
        forall (i=1:(this%nd  )) this%r(i)  = ( rd + ru ) / 2 - cos((2*i-1)*pi/(2*this%nd))/cos(pi/(2*this%nd))/2
        forall (i=1:(this%nd+1)) this%rr(i) = ( rd + ru ) / 2 - cos((  i-1)*pi/(  this%nd))/cos(pi/(2*this%nd))/2
      
      case('homog')
        dr = (ru-rd)/(this%nd-1)
          forall (i=1:(this%nd  )) this%r(i)  = rd        + (i-1)*dr
          forall (i=1:(this%nd+1)) this%rr(i) = rd - dr/2 + (i-1)*dr
    end select

    this%volume = 4 * pi * ( ru**3 - rd**3 ) / 3

  end subroutine init_grid_sub

  subroutine deallocate_grid_sub(this)
    class(T_radialGrid), intent(inout) :: this

    deallocate(this%r, this%rr)

  end subroutine deallocate_grid_sub

end module RadialGrid