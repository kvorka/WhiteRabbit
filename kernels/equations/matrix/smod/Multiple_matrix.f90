submodule(Matrix) Multiple_matrix
  implicit none; contains
  
  module pure complex(kind=dbl) function matrix_multiple_fn(this, i, vector)
    class(T_matrix),   intent(in) :: this
    integer,           intent(in) :: i
    complex(kind=dbl), intent(in) :: vector(:)
    integer                       :: j, k
    
    matrix_multiple_fn = czero
    
    k = i-this%ld-1
      do concurrent ( j = max(1,1-k):min(this%ldu,this%n-k) )
        matrix_multiple_fn = matrix_multiple_fn + this%M(j,i) * vector(j+k)
      end do
    
  end function matrix_multiple_fn
  
end submodule Multiple_matrix