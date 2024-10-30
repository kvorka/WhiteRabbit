submodule (matrix) mult
  implicit none; contains
  
  module procedure matrix_multiple_fn
    integer :: k, indstart, indend
    
    k        = i-this%ld-1
    indstart = max(1,1-k)
    indend   = min(this%ldu,this%n-k)
    
    matrix_multiple_fn = sum( this%M(indstart:indend,i) * vector(k+indstart:k+indend) )
    
  end procedure matrix_multiple_fn
  
end submodule mult