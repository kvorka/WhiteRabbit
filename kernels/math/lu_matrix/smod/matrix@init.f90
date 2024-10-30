submodule (matrix) init
  implicit none; contains
  
  module procedure init_matrix_sub
    
    this%n   = n
    this%ld  = ld
    this%lu  = lu
    this%ldu = ld+1+lu
    
    allocate( this%M( this%ldu, this%n ) ) ; this%M = zero
    allocate( this%U( this%ldu, this%n ) ) ; this%U = zero
    allocate( this%L( this%ld , this%n ) ) ; this%L = zero
    allocate( this%I( this%n )           ) ; this%I = 0
    
  end procedure init_matrix_sub
  
  module procedure deallocate_matrix_sub
    
    deallocate( this%M, this%U, this%L, this%I )
    
  end procedure deallocate_matrix_sub
  
end submodule init