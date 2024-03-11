submodule (Fourier_transform) fft_c2r
  implicit none; contains
  
  module pure subroutine fft_c2r_exec_sub(this, m, cx, x)
    class(T_fft),      intent(in)  :: this
    integer,           intent(in)  :: m
    complex(kind=dbl), intent(in)  :: cx(m,*)
    real(kind=dbl),    intent(out) :: x(m,*)
    integer                        :: i1, i2
    
    do concurrent ( i2 = 1:this%n , i1 = 1:m )
      x(i1,i2) = 0._dbl
    end do
    
    do concurrent (i2 = 1:this%np, i1 = 1:m)
      x(i1,2*i2-1) = cx(i1,i2)%re
      x(i1,2*i2  ) = cx(i1,i2)%im
    end do
    
    call this%fft_c2r_sub( m,  x )
    
  end subroutine fft_c2r_exec_sub
  
  module pure subroutine fft_c2r_sub(this, m, x)
    class(T_fft),   intent(in)    :: this
    integer,        intent(in)    :: m
    real(kind=dbl), intent(inout) :: x(m,2,0:this%n/2-1)
    integer                       :: iv, i, ii, isj, j, isj2
    real(kind=dbl)                :: temp, tempre1, tempim1, tempre2, tempim2
    real(kind=dbl), allocatable   :: y(:,:)
    
    do concurrent ( iv = 1:m )
      temp      = x(iv,1,0)
      x(iv,1,0) = x(iv,1,0) + x(iv,2,0)
      x(iv,2,0) = temp      - x(iv,2,0)
    end do
    
    do concurrent ( i = 1:(this%n-2)/4 , iv = 1:m )
      tempre1 = x(iv,1,i) + x(iv,1,this%n/2-i)
      tempim1 = x(iv,1,i) - x(iv,1,this%n/2-i)
      tempre2 = x(iv,2,i) + x(iv,2,this%n/2-i)
      tempim2 = x(iv,2,i) - x(iv,2,this%n/2-i)
      
      x(iv,1,         i) =       tempre1 - tempim1 * this%t(this%n+2*i  ) - tempre2 * this%t(this%n+2*i-1)
      x(iv,2,         i) =       tempim2 + tempim1 * this%t(this%n+2*i-1) - tempre2 * this%t(this%n+2*i  )
      x(iv,1,this%n/2-i) =   2 * tempre1 - x(iv,1,i)
      x(iv,2,this%n/2-i) = - 2 * tempim2 + x(iv,2,i)
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
        do concurrent ( ii = 1:2, i = 1:m )
          y(i,ii) = x(i,ii,isj)
        end do
        
        do
           j = j + 1 ; isj2 = this%it(j)
           
           if ( isj2 < 0 ) then
              do concurrent ( ii = 1:2, i = 1:m )
                x(i,ii,isj     ) = x(i,ii,isj2-imm)
                x(i,ii,isj2-imm) = y(i,ii)
              end do
              
              j = j + 1 ; exit
           else
              do concurrent ( ii = 1:2, i = 1:m )
                x(i,ii,isj) = x(i,ii,isj2)
              end do
              
              isj = isj2
           end if
        end do
      end if
    end do
    
    deallocate( y )
    
  end subroutine fft_c2r_sub
  
end submodule fft_c2r