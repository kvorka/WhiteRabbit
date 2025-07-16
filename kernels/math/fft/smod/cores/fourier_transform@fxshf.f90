submodule (fourier_transform) fxshf
  implicit none; contains
  
  module procedure fxzshf
    integer                             :: j, i1, i2, i30, i31
    real(kind=dbl), pointer, contiguous :: y(:,:)
    type(c_ptr)                         :: c_y
    
    call alloc_aligned2d_sub( 16, m, c_y, y )
    
    j = 1  
      do while (j <= this%n/2-2)
        i30 = this%it(j)
        
        if (i30 < 0) then
          j = j + 1
          
        else
          do i2 = 1, m
            !$omp simd
            do i1 = 1, 16
              y(i1,i2) = x(i1,i2,i30)
            end do
          end do
          
          do
            j    = j + 1
            i31 = this%it(j)
            
            if ( i31 < 0 ) then
              do i2 = 1, m
                !$omp simd
                do i1 = 1, 16
                  x(i1,i2,i30    ) = x(i1,i2,i31-imm)
                  x(i1,i2,i31-imm) = y(i1,i2)
                end do
              end do
              
              j = j + 1
              exit
            
            else
              do i2 = 1, m
                !$omp simd
                do i1 = 1, 16
                  x(i1,i2,i30) = x(i1,i2,i31)
                end do
              end do
              
              i30 = i31
            end if
          end do
        end if
      end do
    
    call free_aligned2d_sub( c_y, y )
    
  end procedure fxzshf
  
end submodule fxshf