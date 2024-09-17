module RadialGrid
  use Math
  implicit none
  
  type, public :: T_radialGrid
    integer                     :: nd
    real(kind=dbl)              :: volume
    real(kind=dbl), allocatable :: r(:), rr(:)
    
    contains
    
    procedure :: init_sub       => init_grid_sub
    procedure :: deallocate_sub => deallocate_grid_sub
    procedure :: hd, d, c
    procedure :: hdd, dd, cc, hdrr, drr
    procedure :: interpolation_fn
    procedure :: volumetric_integral_real_fn, volumetric_integral_cmplx_fn
    procedure :: radial_integral_real_fn, radial_integral_cmplx_fn
    
    generic :: intV_fn => volumetric_integral_real_fn, volumetric_integral_cmplx_fn
    generic :: intR_fn => radial_integral_real_fn    , radial_integral_cmplx_fn
      
  end type T_radialGrid
  
  interface
    module pure subroutine init_grid_sub(this, nr, rd, ru, grid_type)
      class(T_radialGrid), intent(inout) :: this
      integer,             intent(in)    :: nr
      real(kind=dbl),      intent(in)    :: rd, ru
      character(len=5),    intent(in)    :: grid_type
    end subroutine init_grid_sub
    
    module pure subroutine deallocate_grid_sub(this)
      class(T_radialGrid), intent(inout) :: this
    end subroutine deallocate_grid_sub
    
    module pure real(kind=dbl) function hd(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function hd
    
    module pure real(kind=dbl) function d(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function d
    
    module pure real(kind=dbl) function c(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function c
    
    module pure real(kind=dbl) function hdd(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function hdd
    
    module pure real(kind=dbl) function dd(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function dd
    
    module pure real(kind=dbl) function cc(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function cc
    
    module pure real(kind=dbl) function hdrr(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function hdrr
    
    module pure real(kind=dbl) function drr(this, i, p)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, p
    end function drr
    
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
    
    module pure function interpolation_fn(this, dimOut, i, rr1, field) result(resField)
      class(T_radialGrid), intent(in) :: this
      integer,             intent(in) :: i, dimOut
      real(kind=dbl),      intent(in) :: rr1(:)
      complex(kind=dbl),   intent(in) :: field(:,:)
      complex(kind=dbl), allocatable  :: resField(:)
    end function interpolation_fn
  end interface
  
end module RadialGrid