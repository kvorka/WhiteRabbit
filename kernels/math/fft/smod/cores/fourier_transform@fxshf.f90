submodule (fourier_transform) fxshf
  implicit none; contains
  
  module procedure fxzshf
    integer                             :: i1, iv, j, isj, isj2
    real(kind=dbl), pointer, contiguous :: y(:,:)
    type(c_ptr)                         :: c_y
    
    call alloc_aligned2d_sub( 16, m, c_y, y )
    
    j = 1  
      do while (j <= this%n/2-2)
        isj = this%it(j)
        
        if (isj < 0) then
          j = j + 1
          
        else
          do iv = 1, m
            !$omp simd
            do i1 = 1, 16
              y(i1,iv) = x(i1,iv,isj)
            end do
          end do
          
          do
            j    = j + 1
            isj2 = this%it(j)
            
            if ( isj2 < 0 ) then
              do iv = 1, m
                !$omp simd
                do i1 = 1, 16
                  x(i1,iv,isj)      = x(i1,iv,isj2-imm)
                  x(i1,iv,isj2-imm) = y(i1,iv)
                end do
              end do
              
              j = j + 1
              exit
            
            else
              do iv = 1, m
                !$omp simd
                do i1 = 1, 16
                  x(i1,iv,isj) = x(i1,iv,isj2)
                end do
              end do
              
              isj = isj2
            end if
          end do
        end if
      end do
    
    call free_aligned2d_sub( c_y, y )
    
  end procedure fxzshf
  
end submodule fxshf