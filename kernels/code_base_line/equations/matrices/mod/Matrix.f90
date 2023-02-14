module Matrix
  use Math
  use LU
  implicit none
  
  type, public :: T_matrix
    
    integer                     :: n, ns, nu               !dimenzia matice, rozmer pod/naddiagonaly
    real(kind=dbl), allocatable :: M(:,:), U(:,:), L(:,:)  !povodna, horna a dolna trojuholnikova matica sustavy
    integer,        allocatable :: I(:)                    !zaznam o pivotacii pri LU rozklade

    contains

    procedure :: init_sub       => initMatrix_sub
    procedure :: fill_sub       => fillMatrix_sub
    procedure :: luSolve_sub    => luSolutionMatrix_sub
    procedure :: multipl_fn     => matrixMultiple_fn
    procedure :: multipl2_fn    => matrixMultiple2_fn
    procedure :: deallocate_sub => deallocate_sub

  end type T_matrix
  
  private :: initMatrix_sub
  private :: fillMatrix_sub
  private :: luSolutionMatrix_sub
  private :: matrixMultiple_fn
  private :: matrixMultiple2_fn
  private :: deallocate_sub

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

  pure function matrixMultiple2_fn(this, i, vector) result(multipl)
    class(T_matrix),   intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl), intent(in) :: vector(:,:)
    complex(kind=dbl)             :: multipl(size(vector,2))

    multipl = multipl2_fn(i, this%n, this%ns, this%nu, this%M, vector)

  end function matrixMultiple2_fn

  subroutine deallocate_sub(this)
    class(T_matrix), intent(inout) :: this

    deallocate(this%M, this%U, this%L, this%I)

  end subroutine deallocate_sub

end module Matrix
