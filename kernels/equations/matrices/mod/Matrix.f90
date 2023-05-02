module Matrix
  use LU
  implicit none
  
  type, public :: T_matrix
    integer                     :: n, ns, nu
    real(kind=dbl), allocatable :: M(:,:), U(:,:), L(:,:)
    integer,        allocatable :: I(:)

    contains

    procedure, public, pass :: init_sub       => initMatrix_sub
    procedure, public, pass :: fill_sub       => fillMatrix_sub
    procedure, public, pass :: luSolve_sub    => luSolutionMatrix_sub
    procedure, public, pass :: multipl_fn     => matrixMultiple_fn
    procedure, public, pass :: deallocate_sub => deallocateMatrix_sub

  end type T_matrix

  contains
  
  subroutine initMatrix_sub(this, n, ns, nu)
    class(T_matrix), intent(inout) :: this
    integer,         intent(in)    :: n, ns, nu

    allocate(this%M(ns+1+nu, n), this%U(ns+1+nu, n), this%L(ns, n), this%I(n))

    this%n = n; this%ns = ns; this%nu = nu
    this%M = 0._dbl; this%U = 0._dbl; this%L = 0._dbl; this%I = 0

  end subroutine initMatrix_sub

  subroutine fillMatrix_sub(this, matrixU, matrixM)
    class(T_matrix), intent(inout) :: this
    real(kind=dbl),  intent(in)    :: matrixU(:,:), matrixM(:,:)

    this%M = matrixM; this%U = matrixU
      call ludecomposition_sub(this%n, this%ns, this%nu, this%U, this%L, this%I)

  end subroutine fillMatrix_sub
  
  subroutine luSolutionMatrix_sub(this, RHS)
    class(T_matrix),   intent(in)    :: this
    complex(kind=dbl), intent(inout) :: RHS(:)

    call lusolution_sub(this%n, this%ns, this%nu, this%U, this%L, this%I, RHS)

  end subroutine luSolutionMatrix_sub

  pure complex(kind=dbl) function matrixMultiple_fn(this, i, vector)
    class(T_matrix),   intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl), intent(in) :: vector(:)

    matrixMultiple_fn = multipl1_fn(i, this%n, this%ns, this%nu, this%M, vector)

  end function matrixMultiple_fn

  subroutine deallocateMatrix_sub(this)
    class(T_matrix), intent(inout) :: this

    deallocate(this%M, this%U, this%L, this%I)

  end subroutine deallocateMatrix_sub

end module Matrix