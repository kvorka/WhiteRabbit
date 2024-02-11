submodule(Matrix) Init_matrix
  implicit none; contains
  
  module pure subroutine init_matrix_sub(this, n, ld, lu)
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
    
  end subroutine init_matrix_sub
  
  module pure subroutine deallocate_matrix_sub(this)
    class(T_matrix), intent(inout) :: this
    
    deallocate( this%M, this%U, this%L, this%I )
    
  end subroutine deallocate_matrix_sub
  
end submodule Init_matrix