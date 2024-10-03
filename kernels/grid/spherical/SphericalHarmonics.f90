module SphericalHarmonics
  use nulify
  use legendre_function
  use legendre_polynom
  use fourier_transform
  use sphsvt
  implicit none
  
  integer, parameter :: addmissible_jmax(48) = [   5,   7,   9,  13,  15,  21,  27,  29,  33,  37,  45, 47,  51,  57,  61,  69,  &
                                               &  77,  87,  93,  97, 105, 117, 125, 141, 147, 157, 159, 177, 189, 197, 213, 237, &
                                               & 247, 253, 267, 285, 297, 317, 321, 357, 381, 397, 429, 447, 477, 483, 497, 997  ]
  
  type, public :: T_lateralGrid
    integer,                     public  :: jmax, nLegendre, nFourier
    real(kind=dbl), allocatable, public  :: cosx(:), weight(:), amj(:), bmj(:), cmm(:)
    type(T_fft),                 public  :: fourtrans
    type(T_sphsvt),              public  :: reindexing
    
    contains
    
    procedure :: init_sub       => init_harmonics_sub
    procedure :: deallocate_sub => deallocate_harmonics_sub
    
    procedure :: transform_sub
    procedure :: vcss_sub, vcsv_sub, vcst_sub, vcvv_sub, vcvgv_sub, vcvv_vcvgv_sub, vcss_add_vcvv_sub
    procedure :: grid_to_space_sub, space_to_grid_sub
    
  end type T_lateralGrid
  
  interface
    module subroutine init_harmonics_sub(this, jmax)
      class(T_lateralGrid), intent(inout) :: this
      integer,              intent(in)    :: jmax
    end subroutine init_harmonics_sub
    
    module pure subroutine deallocate_harmonics_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_harmonics_sub
    
    module pure subroutine transform_sub(this, nf, nb, cc, cr, grid_sub)
      class(T_lateralGrid), intent(in)    :: this
      integer,              intent(in)    :: nf, nb
      complex(kind=dbl),    intent(in)    :: cc(nb,*)
      complex(kind=dbl),    intent(inout) :: cr(nf,*)
      
      interface
        module pure subroutine grid_sub(step, nfour, gxyz)
          integer,                intent(in)    :: nfour, step
          real(kind=dbl), target, intent(inout) :: gxyz(*)
        end subroutine grid_sub
      end interface
    end subroutine transform_sub
    
    module pure subroutine space_to_grid_sub(this, cc, grid)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cc(*)
      real(kind=dbl),       intent(out) :: grid(:,:,:)
    end subroutine space_to_grid_sub
    
    module pure subroutine grid_to_space_sub(this, grid, cr)
      class(T_lateralGrid), intent(in)    :: this
      real(kind=dbl),       intent(inout) :: grid(:,:,:)
      complex(kind=dbl),    intent(out)   :: cr(*)
    end subroutine grid_to_space_sub
    
    module pure subroutine vcss_sub(this, cajm, cbjm, cjm)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcss_sub
    
    module pure subroutine vcsv_sub(this, cajm, cbjml, cjm)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cajm(*), cbjml(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcsv_sub
    
    module pure subroutine vcst_sub(this, cajm, cbjml2, cjml2)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cajm(*), cbjml2(*)
      complex(kind=dbl),    intent(out) :: cjml2(*)
    end subroutine vcst_sub
    
    module pure subroutine vcvv_sub(this, cajml, cbjml, cjm)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cajml(*), cbjml(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcvv_sub
    
    module pure subroutine vcvgv_sub(this, ri, dv_r, v, cjm)
      class(T_lateralGrid), intent(in)  :: this
      real(kind=dbl),       intent(in)  :: ri
      complex(kind=dbl),    intent(in)  :: v(*), dv_r(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcvgv_sub
    
    module pure subroutine vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
      class(T_lateralGrid), intent(in)  :: this
      real(kind=dbl),       intent(in)  :: ri
      complex(kind=dbl),    intent(in)  :: dv_r(*), q(*), v(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcvv_vcvgv_sub
    
    module pure subroutine vcss_add_vcvv_sub(this, cajm, cbjm, cajml, cbjml, cjm)
      class(T_lateralGrid), intent(in)  :: this
      complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*), cajml(*), cbjml(*)
      complex(kind=dbl),    intent(out) :: cjm(*)
    end subroutine vcss_add_vcvv_sub
  end interface
  
end module SphericalHarmonics
