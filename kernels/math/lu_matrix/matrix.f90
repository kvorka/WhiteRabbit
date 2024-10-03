module matrix
  !Original work on LU: Numerical Recipes in Fortran77
  use math
  implicit none
  
  type, public :: T_matrix
    integer                     :: n, ld, lu, ldu
    real(kind=dbl), allocatable :: M(:,:), U(:,:), L(:,:)
    integer,        allocatable :: I(:)
    
    contains
    
    procedure :: init_sub       => init_matrix_sub
    procedure :: fill_sub       => lu_decomposition_sub
    procedure :: luSolve_sub    => lu_solve_sub
    procedure :: multipl_fn     => matrix_multiple_fn
    procedure :: deallocate_sub => deallocate_matrix_sub
    
  end type T_matrix
  
  interface
    module pure subroutine init_matrix_sub(this, n, ld, lu)
      class(T_matrix), intent(inout) :: this
      integer,         intent(in)    :: n, ld, lu
    end subroutine init_matrix_sub
    
    module pure subroutine deallocate_matrix_sub(this)
      class(T_matrix), intent(inout) :: this
    end subroutine deallocate_matrix_sub
    
    module pure subroutine lu_decomposition_sub(this, matrixU, matrixM)
      class(T_matrix), intent(inout) :: this
      real(kind=dbl),  intent(in)    :: matrixU(:,:), matrixM(:,:)
    end subroutine lu_decomposition_sub
    
    module pure subroutine lu_solve_sub(this, b)
      class(T_matrix),   intent(in)    :: this
      complex(kind=dbl), intent(inout) :: b(:)
    end subroutine lu_solve_sub
    
    module pure complex(kind=dbl) function matrix_multiple_fn(this, i, vector)
      class(T_matrix),   intent(in) :: this
      integer,           intent(in) :: i
      complex(kind=dbl), intent(in) :: vector(:)
    end function matrix_multiple_fn
  end interface
  
end module matrix