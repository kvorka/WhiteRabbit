module Matrix
  !Original work on LU: Numerical Recipes in Fortran77
  use Math
  implicit none
  
  type, public :: T_matrix
    integer                     :: n, ld, lu, ldu
    real(kind=dbl), allocatable :: M(:,:), U(:,:), L(:,:)
    integer,        allocatable :: I(:)
    
    contains
    
    procedure :: init_sub       => initMatrix_sub
    procedure :: fill_sub       => luDecompMatrix_sub
    procedure :: luSolve_sub    => luSolutionMatrix_sub
    procedure :: multipl_fn     => matrixMultiple_fn
    procedure :: deallocate_sub => deallocateMatrix_sub
    
  end type T_matrix
  
  contains
  
  pure subroutine initMatrix_sub(this, n, ld, lu)
    class(T_matrix), intent(inout) :: this
    integer,         intent(in)    :: n, ld, lu
    
    this%n   = n
    this%ld  = ld
    this%lu  = lu
    this%ldu = ld+1+lu
    
    allocate( this%M( this%ldu, this%n ) ) ; this%M = zero
    allocate( this%U( this%ldu, this%n ) ) ; this%U = zero
    allocate( this%L( this%ld , this%n ) ) ; this%L = zero
    allocate( this%I( this%n )           ) ; this%I = 0
    
  end subroutine initMatrix_sub
  
  pure subroutine luDecompMatrix_sub(this, matrixU, matrixM)
    class(T_matrix), intent(inout) :: this
    real(kind=dbl),  intent(in)    :: matrixU(:,:), matrixM(:,:)
    integer                        :: i, j, k
    real(kind=dbl), allocatable    :: pom(:)
    
    this%M = matrixM
    this%U = matrixU
    
    k = this%ld
      do i = 1, this%ld
        this%U(this%ld+2-i-k:this%ldu-k,i) = this%U(this%ld+2-i:this%ldu,i)
        
        k = k-1
          this%U(this%ldu-k:this%ldu,i) = zero
      end do
    
    allocate( pom(this%ldu) )
    
    do j = 1, this%n
      k         = min(this%ld+j,this%n)
      i         = maxloc( abs( this%U(1,j:k) ), 1 ) + j - 1
      this%I(j) = i
      
      if (i /= j) then
        pom         = this%U(:,j)
        this%U(:,j) = this%U(:,i)
        this%U(:,i) = pom
      end if
      
      pom = this%U(:,j)
        do i = j+1, k
          this%L(i-j,j)          = this%U(1,i) / pom(1)
          this%U(1:this%ldu-1,i) = this%U(2:this%ldu,i) - this%L(i-j,j) * pom(2:this%ldu)
          this%U(this%ldu,i)     = zero
        end do
    end do
    
    deallocate(pom)
    
  end subroutine luDecompMatrix_sub
  
  pure subroutine luSolutionMatrix_sub(this, b)
    class(T_matrix),   intent(in)    :: this
    complex(kind=dbl), intent(inout) :: b(:)
    integer                          :: i, j, k
    complex(kind=dbl)                :: dum
    
    do j = 1, this%n
      i = this%I(j) ; dum = b(i)
        if (i /= j) then
          b(i) = b(j)
          b(j) = dum
        end if
      
      k = min(this%n,this%ld+j) ; b(j+1:k) = b(j+1:k) - this%L(1:k-j,j) * dum
    end do
    
    do i = this%n, 1, -1
      k = min(this%ldu,this%n-i+1) ; b(i) = ( b(i) - sum( this%U(2:k,i) * b(i+1:i+k-1) ) ) / this%U(1,i)
    end do
    
  end subroutine luSolutionMatrix_sub
  
  pure complex(kind=dbl) function matrixMultiple_fn(this, i, vector)
    class(T_matrix),   intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl), intent(in) :: vector(:)
    integer                       :: indM1, indM2, indV1, indV2
    
    indM1 = max( 1 ,   this%ld+2-i ) ; indM2 = min( this%ldu  , this%ld+1+this%n-i )
    indV1 = max( 1 , i-this%ld     ) ; indV2 = min( i+this%lu ,           this%n   )
    
    matrixMultiple_fn = sum( this%M( indM1 : indM2, i ) * vector( indV1 : indV2 ) )
    
  end function matrixMultiple_fn
  
  pure subroutine deallocateMatrix_sub(this)
    class(T_matrix), intent(inout) :: this
    
    deallocate( this%M, this%U, this%L, this%I )
    
  end subroutine deallocateMatrix_sub
  
end module Matrix