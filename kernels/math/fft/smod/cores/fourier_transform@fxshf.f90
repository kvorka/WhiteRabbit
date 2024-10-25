submodule (fourier_transform) fxshf
  implicit none; contains
  
  module pure subroutine fxzshf(this, m, x)
    class(T_fft),   intent(in)    :: this
    integer,        intent(in)    :: m
    real(kind=dbl), intent(inout) :: x(m,2,0:this%n/2-1)
    integer                       :: i, iv, j, isj, isj2
    real(kind=dbl), allocatable   :: y(:,:)
    
    allocate( y(m,2) )
    
    j = 1  
      do while (j <= this%n/2-2)
        isj = this%it(j)
        
        if (isj < 0) then
          j = j + 1
        else
          do concurrent ( i = 1:2, iv = 1:m )
            y(iv,i) = x(iv,i,isj)
          end do
          
          do
            j    = j + 1
            isj2 = this%it(j)
            
            if ( isj2 < 0 ) then
              do concurrent ( i = 1:2, iv = 1:m )
                x(iv,i,isj)      = x(iv,i,isj2-imm)
                x(iv,i,isj2-imm) = y(iv,i)
              end do
              
              j = j + 1
              exit
            else
              do concurrent ( i = 1:2, iv = 1:m )
                x(iv,i,isj) = x(iv,i,isj2)
              end do
              
              isj = isj2
            end if
          end do
        end if
      end do
    
    deallocate( y )
    
  end subroutine fxzshf
  
end submodule fxshf