submodule(Matrix) LU_solve
  implicit none; contains
  
  module pure subroutine lu_solve_sub(this, b)
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
      
      do concurrent ( i = j+1:min(this%n,this%ld+j) )
        b(i) = b(i) - this%L(i-j,j) * dum
      end do
    end do
    
    do i = this%n, 1, -1
      k = min(this%ldu,this%n-i+1)
        b(i) = ( b(i) - sum( this%U(2:k,i) * b(i+1:i+k-1) ) ) / this%U(1,i)
    end do
    
  end subroutine lu_solve_sub
  
end submodule LU_solve