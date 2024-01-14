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
    integer                            :: i, j, k, l
    real(kind=dbl)                     :: pom
    
    k = this%ld
      do i = 1, this%ld
        do concurrent ( l = this%ld+2-i-k:this%ldu-k )
          this%U(l,i) = this%U(l+k,i)
        end do
        
        k = k-1
          do concurrent ( l = this%ldu-k:this%ldu )
            this%U(l,i) = 0._dbl
          end do
      end do
    
    do j = 1, this%n
      k = min(this%ld+j,this%n)
      
      i = j; pom = abs( this%U(1,j) )
        do l = j, k
          if ( abs(this%U(1,l)) > pom ) then
            i   = l
            pom = abs(this%U(1,l))
          end if
        end do
      
      this%I(j) = i
        if (i /= j) then
          do concurrent ( l = 1:this%ldu )
            pom         = this%U(l,j)
            this%U(l,j) = this%U(l,i)
            this%U(l,i) = pom
          end do
        end if
      
      do i = j+1, k
        this%L(i-j,j) = this%U(1,i) / this%U(1,j)
        pom           = this%L(i-j,j)
        
        do l = 1, this%ldu-1
          this%U(l,i) = this%U(l+1,i) - pom * this%U(l+1,j)
        end do
        
        this%U(this%ldu,i) = 0._dbl
      end do
    end do
    
  end subroutine luDecompMatrix_sub
  
  pure subroutine luSolutionMatrix_sub(this, b)
    class(T_matrix),   intent(in)    :: this
    complex(kind=dbl), intent(inout) :: b(:)
    integer                            :: i, j, k
    complex(kind=dbl)                  :: dum
    
    do j = 1, this%n
      i = this%I(j) ; dum = b(i)
        if (i /= j) then
          b(i) = b(j)
          b(j) = dum
        end if
      
      do concurrent ( i = j+1:min(this%n,this%ld+j) )
        b(i) = b(i) - this%L(i-j,j) * dum
      end do
    end do
    
    do i = this%n, 1, -1
      k = min(this%ldu,this%n-i+1)
        b(i) = ( b(i) - sum( this%U(2:k,i) * b(i+1:i+k-1) ) ) / this%U(1,i)
    end do
    
  end subroutine luSolutionMatrix_sub
  
  pure complex(kind=dbl) function matrixMultiple_fn(this, i, vector)
    class(T_matrix),   intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl), intent(in) :: vector(:)
    integer                       :: j, k
    
    matrixMultiple_fn = czero
    
    k = i-this%ld-1
      do concurrent ( j = max(1,1-k):min(this%ldu,this%n-k) )
        matrixMultiple_fn = matrixMultiple_fn + this%M(j,i) * vector(j+k)
      end do
    
  end function matrixMultiple_fn
  
  pure subroutine deallocateMatrix_sub(this)
    class(T_matrix), intent(inout) :: this
    
    deallocate( this%M, this%U, this%L, this%I )
    
  end subroutine deallocateMatrix_sub
  
end module Matrix