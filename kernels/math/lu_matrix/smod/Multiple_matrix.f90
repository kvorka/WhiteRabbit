submodule(Matrix) Multiple_matrix
  implicit none; contains
  
  module pure complex(kind=dbl) function matrix_multiple_fn(this, i, vector)
    class(T_matrix),   intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl), intent(in) :: vector(:)
    integer                       :: k, indstart, indend
    
    k        = i-this%ld-1
    indstart = max(1,1-k)
    indend   = min(this%ldu,this%n-k)
    
    matrix_multiple_fn = sum( this%M(indstart:indend,i) * vector(k+indstart:k+indend) )
    
  end function matrix_multiple_fn
  
end submodule Multiple_matrix