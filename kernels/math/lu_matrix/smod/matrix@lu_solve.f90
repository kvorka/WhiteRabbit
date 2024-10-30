submodule (matrix) lu_solve
  implicit none; contains
  
  module procedure lu_solve_sub
    integer           :: i, j
    complex(kind=dbl) :: dum
    
    do j = 1, this%n
      i   = this%I(j)
      dum = b(i)
      
      if (i /= j) then
        b(i) = b(j)
        b(j) = dum
      end if
      
      do concurrent ( i = j+1:min(this%n,this%ld+j) )
        b(i) = b(i) - this%L(i-j,j) * dum
      end do
    end do
    
    do i = this%n, 1, -1
      dum = b(i)
        do j = 2, min(this%ldu,this%n-i+1)
          dum = dum - this%U(j,i) * b(i+j-1)
        end do
      
      b(i) = dum / this%U(1,i)
    end do
    
  end procedure lu_solve_sub
  
end submodule lu_solve
