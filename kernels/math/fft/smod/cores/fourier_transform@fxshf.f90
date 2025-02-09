submodule (fourier_transform) fxshf
  implicit none; contains
  
  module procedure fxzshf
    integer                     :: i, iv, j, isj, isj2
    real(kind=dbl), allocatable :: y(:,:)
    
    allocate( y(m,2) )
    
    j = 1  
      do while (j <= this%n/2-2)
        isj = this%it(j)
        
        if (isj < 0) then
          j = j + 1
        else
          !$omp simd collapse(2)
          do i = 1, 2
            do iv = 1, m
              y(iv,i) = x(iv,i,isj)
            end do
          end do
          
          do
            j    = j + 1
            isj2 = this%it(j)
            
            if ( isj2 < 0 ) then
              !$omp simd collapse(2)
              do i = 1, 2
                do iv = 1, m
                  x(iv,i,isj)      = x(iv,i,isj2-imm)
                  x(iv,i,isj2-imm) = y(iv,i)
                end do
              end do
              
              j = j + 1
              exit
            else
              !$omp simd collapse(2)
              do i = 1, 2
                do iv = 1, m
                  x(iv,i,isj) = x(iv,i,isj2)
                end do
              end do
              
              isj = isj2
            end if
          end do
        end if
      end do
    
    deallocate( y )
    
  end procedure fxzshf
  
end submodule fxshf