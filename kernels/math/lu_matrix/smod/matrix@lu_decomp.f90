submodule (matrix) lu_decomp
  implicit none; contains
  
  module procedure lu_decomposition_sub
    integer        :: i, j, k, l
    real(kind=dbl) :: pom
    
    this%U = matrixU
    this%M = matrixM
    
    k = this%ld
      do i = 1, this%ld
        do concurrent ( l = this%ld+2-i-k:this%ldu-k )
          this%U(l,i) = this%U(l+k,i)
        end do
        
        k = k-1
          do concurrent ( l = this%ldu-k:this%ldu )
            this%U(l,i) = zero
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
        
        this%U(this%ldu,i) = zero
      end do
    end do
    
  end procedure lu_decomposition_sub
  
end submodule lu_decomp