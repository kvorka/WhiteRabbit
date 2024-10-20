submodule (fourier_transform) c2r
  implicit none ; contains
  
  module pure subroutine fft_c2r_sub(this, m, x)
    class(T_fft),   intent(in)    :: this
    integer,        intent(in)    :: m
    real(kind=dbl), intent(inout) :: x(m,2,0:this%n/2-1)
    integer                       :: iv, i, isj, j, isj2
    real(kind=dbl)                :: tempre1, tempim1, tempre2, tempim2, t1, t2
    real(kind=dbl), allocatable   :: y(:,:)
    
    do concurrent ( iv = 1:m )
      tempre1   = x(iv,1,0)
      x(iv,1,0) = x(iv,1,0) + x(iv,2,0)
      x(iv,2,0) = tempre1   - x(iv,2,0)
    end do
    
    do i = 1, (this%n-2)/4
      t1 = this%t(this%n+2*i-1)
      t2 = this%t(this%n+2*i  )
      
      do concurrent ( iv = 1:m )
        tempre1 = x(iv,1,i) + x(iv,1,this%n/2-i)
        tempim1 = x(iv,1,i) - x(iv,1,this%n/2-i)
        tempre2 = x(iv,2,i) + x(iv,2,this%n/2-i)
        tempim2 = x(iv,2,i) - x(iv,2,this%n/2-i)
        
        x(iv,1,         i) =       tempre1 - tempim1 * t2 - tempre2 * t1
        x(iv,2,         i) =       tempim2 + tempim1 * t1 - tempre2 * t2
        x(iv,1,this%n/2-i) =   2 * tempre1 - x(iv,1,i)
        x(iv,2,this%n/2-i) = - 2 * tempim2 + x(iv,2,i)
      end do
    end do
    
    if ( mod(this%n,4) == 0 ) then
      do concurrent ( iv = 1:m )
        x(iv,1,this%n/4) = +2 * x(iv,1,this%n/4)
        x(iv,2,this%n/4) = -2 * x(iv,2,this%n/4)
      end do
    end if
    
    select case ( mod(this%it(this%n/2-1),4)+2 )
      case (4)
        call fxzm4b(m, this%n/2, x)
      case (2)
        call fxzm2b(m, this%n/2, x)
      case (3)
        call fxzm3b(m, this%n/2, x)
      case (5)
        call fxzm5b(m, this%n/2, x)
    end select
    
    call fxztal(m, 1, this%n/2, x, this%t, 1, this%it(this%n/2), 0, 0, this%it(this%n/2-1))
    
    allocate( y(m,2) ) ; j = 1
      
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
            
            j = j + 1 ; exit
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
    
  end subroutine fft_c2r_sub
  
end submodule c2r