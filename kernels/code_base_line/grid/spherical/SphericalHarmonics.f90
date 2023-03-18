module SphericalHarmonics
  use Clebsch_Legendre
  use, intrinsic :: iso_c_binding
  implicit none
  include '/usr/include/fftw3.f03'

  integer, parameter :: fftw_flags = FFTW_DESTROY_INPUT + FFTW_EXHAUSTIVE
  integer, parameter :: step       = 8

  type, public :: T_lateralGrid
    integer,                     private :: jmax, jms, jms1, jms2, jmv, jmv1, maxj, nLegendre, nFourier
    real(kind=dbl), allocatable, private :: roots(:), fftLege(:), amjrr(:), bmjrr(:), cmmrr(:)
    type(C_ptr),                 private :: fftw_01_back, fftw_03_back, fftw_04_back
    type(C_ptr),                 private :: fftw_04_forw, fftw_06_forw, fftw_12_forw, fftw_16_forw, fftw_19_forw
  
    contains
  
    procedure :: init_sub       => init_harmonics_sub
    procedure :: deallocate_sub => deallocate_harmonics_sub

    procedure :: init_vcsv_sub
    procedure :: vcsv_fn
    procedure :: deallocate_fftw_vcsv_sub

    procedure :: init_vcvv_sub
    procedure :: vcvv_fn
    procedure :: deallocate_fftw_vcvv_sub

    procedure :: init_vcvgv_sub
    procedure :: vcvgv_fn
    procedure :: deallocate_fftw_vcvgv_sub

    procedure :: init_vcsv_vcvgv_sub
    procedure :: vcsv_vcvgv_fn
    procedure :: deallocate_fftw_vcsv_vcvgv_sub

    procedure :: init_vcsv_vcvv_vcvgv_sub
    procedure :: vcsv_vcvv_vcvgv_sub
    procedure :: deallocate_fftw_vcsv_vcvv_vcvgv_sub

  end type T_lateralGrid

  interface
    module subroutine init_vcsv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine init_vcsv_sub

    module function vcsv_fn(this, cajm, cbjml) result(cjml)
      class(T_lateralGrid), intent(in) :: this
      complex(kind=dbl),    intent(in) :: cajm(:), cbjml(:)
      complex(kind=dbl)                :: cjml(this%jmv)
    end function vcsv_fn

    module subroutine deallocate_fftw_vcsv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_fftw_vcsv_sub
  end interface

  interface
    module subroutine init_vcvv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine init_vcvv_sub

    module function vcvv_fn(this, cajml, cbjml) result(cjm)
      class(T_lateralGrid), intent(in) :: this
      complex(kind=dbl),    intent(in) :: cajml(:), cbjml(:)
      complex(kind=dbl)                :: cjm(this%jmax*(this%jmax+1)/2+this%jmax+1)
    end function vcvv_fn

    module subroutine deallocate_fftw_vcvv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_fftw_vcvv_sub
  end interface

  interface
    module subroutine init_vcvgv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine init_vcvgv_sub

    module function vcvgv_fn(this, cajml, cbjml) result(cjml)
      class(T_lateralGrid), intent(in) :: this
      complex(kind=dbl),    intent(in) :: cajml(:), cbjml(:)
      complex(kind=dbl)                :: cjml(this%jmv)
    end function vcvgv_fn

    module subroutine deallocate_fftw_vcvgv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_fftw_vcvgv_sub
  end interface

  interface
    module subroutine init_vcsv_vcvgv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine init_vcsv_vcvgv_sub

    module function vcsv_vcvgv_fn(this, ri, dv_r, v) result(cjml)
      class(T_lateralGrid), intent(in) :: this
      real(kind=dbl),       intent(in) :: ri
      complex(kind=dbl),    intent(in) :: v(:), dv_r(:)
      complex(kind=dbl)                :: cjml(this%jmv)
    end function vcsv_vcvgv_fn

    module subroutine deallocate_fftw_vcsv_vcvgv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_fftw_vcsv_vcvgv_sub
  end interface

  interface
    module subroutine init_vcsv_vcvv_vcvgv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine init_vcsv_vcvv_vcvgv_sub

    module subroutine vcsv_vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm, cjml)
      class(T_lateralGrid), intent(in)  :: this
      real(kind=dbl),       intent(in)  :: ri
      complex(kind=dbl),    intent(in)  :: dv_r(:), q(:), v(:)
      complex(kind=dbl),    intent(out) :: cjm(:), cjml(:)
    end subroutine vcsv_vcvv_vcvgv_sub
    
    module subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub(this)
      class(T_lateralGrid), intent(inout) :: this
    end subroutine deallocate_fftw_vcsv_vcvv_vcvgv_sub
  end interface

  private :: init_harmonics_sub
  private :: deallocate_harmonics_sub
  private :: destroy_plan_sub

  contains

  subroutine init_harmonics_sub(this, jmax)
    class(T_lateralGrid), intent(inout) :: this
    integer,              intent(in)    :: jmax
    integer                             :: i, k, j, m, n, ncnt
    real(kind=dbl)                      :: xincr, x, y, fx, fy
    integer,           allocatable      :: iemb(:), oemb(:)
    real(kind=dbl),    allocatable      :: testField_re(:,:), testField3_re(:,:,:)
    complex(kind=dbl), allocatable      :: testField(:,:), testField3(:,:,:)

    this%jmax = jmax
    this%jms  =     ( jmax   *(jmax+1)/2 +  jmax   ) + 1
    this%jms1 =     ((jmax+1)*(jmax+2)/2 + (jmax+1)) + 1
    this%jms2 =     ((jmax+2)*(jmax+3)/2 + (jmax+2)) + 1
    this%jmv  = 3 * ( jmax   *(jmax+1)/2 +  jmax   ) + 1
    this%jmv1 = 3 * ((jmax+1)*(jmax+2)/2 + (jmax+1)) + 1
    
    this%maxj      = jmax+2
    this%nFourier  = 3*(this%maxj+1)
    this%nLegendre = (((3*(this%maxj+1)/2+1)/2+1+step)/step)*step

    allocate( this%amjrr(this%jms2), this%bmjrr(this%jms2), this%cmmrr(0:this%maxj), &
            & this%roots(this%nLegendre), this%fftLege(this%nLegendre)               )

    n = this%nLegendre
      do
        n = 2*n; xincr = 1._dbl/n
        if (xincr < 1.0d-15) exit

        ncnt = 0
        x = 0._dbl; fx = lege_fn(2*this%nLegendre, x)
        do i = 1, n
          y = x + xincr; fy = lege_fn(2*this%nLegendre, y)
          if (fx*fy < 0._dbl) ncnt = ncnt+1
          x = y; fx = fy
        end do

        if (ncnt == this%nLegendre) exit
      end do

    i = 0
    x = 0._dbl   ; fx = lege_fn(2*this%nLegendre, x)
    y = x + xincr; fy = lege_fn(2*this%nLegendre, y)
      do
        if (fx*fy < 0._dbl) then
          i = i+1
            this%roots(i) = xnode_fn(this%nLegendre, x, y, fx, fy)
            if (i == this%nLegendre) exit
        end if

        x = y;         fx = fy
        y = x + xincr; fy = lege_fn(2*this%nLegendre, y)
      end do

    do m = 0, this%maxj
      this%cmmrr(m) = sqrt( (2*m+3._dbl) / (m+1._dbl) / 2 )

      do j = m+1, this%maxj
        this%amjrr(m*(this%maxj+1)-m*(m+1)/2+j+1) = sqrt((2*j-1._dbl)*(2*j+1._dbl)                          /(        (j-m)*(j+m)))
        this%bmjrr(m*(this%maxj+1)-m*(m+1)/2+j+1) = sqrt(             (2*j+1._dbl)*(j-m-1._dbl)*(j+m-1._dbl)/((2*j-3)*(j-m)*(j+m)))
      end do
    end do
      
    do i = 1, this%nLegendre
      this%fftLege(i) = (1-this%roots(i)**2) / lege_fn(2*this%nLegendre-1, this%roots(i))**2
    end do
    
  end subroutine init_harmonics_sub

  subroutine deallocate_harmonics_sub(this)
    class(T_lateralGrid), intent(inout) :: this

    call destroy_plan_sub(this%fftw_01_back)
    call destroy_plan_sub(this%fftw_03_back)
    call destroy_plan_sub(this%fftw_04_back)

    call destroy_plan_sub(this%fftw_04_forw)
    call destroy_plan_sub(this%fftw_06_forw)
    call destroy_plan_sub(this%fftw_12_forw)
    call destroy_plan_sub(this%fftw_16_forw)
    call destroy_plan_sub(this%fftw_19_forw)

    deallocate(this%roots, this%amjrr, this%bmjrr, this%cmmrr, this%fftLege)

  end subroutine deallocate_harmonics_sub

  subroutine destroy_plan_sub(plan)
    type(C_ptr), intent(inout) :: plan

    if ( c_associated(plan) ) call fftw_destroy_plan(plan)

  end subroutine destroy_plan_sub
  
end module SphericalHarmonics
