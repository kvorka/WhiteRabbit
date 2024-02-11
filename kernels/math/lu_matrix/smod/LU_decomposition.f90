submodule(Matrix) LU_decomposition
  implicit none; contains
  
  module pure subroutine lu_decomposition_sub(this, matrixU, matrixM)
    class(T_matrix), intent(inout) :: this
    real(kind=dbl),  intent(in)    :: matrixU(:,:), matrixM(:,:)
    integer                        :: i, j, k, l
    real(kind=dbl)                 :: pom
    
    do concurrent ( j = 1:this%n , i = 1:this%ldu )
      this%U(i,j) = matrixU(i,j)
      this%M(i,j) = matrixM(i,j)
    end do
    
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
    
  end subroutine lu_decomposition_sub
  
end submodule LU_decomposition