submodule (Fourier_transform) fft_r2c
  implicit none; contains
  
  module pure subroutine fft_r2c_exec_sub(this, m, x, cx)
    class(T_fft),      intent(in)    :: this
    integer,           intent(in)    :: m
    real(kind=dbl),    intent(inout) :: x(m,2,*)
    complex(kind=dbl), intent(out)   :: cx(m,*)
    integer                          :: i1, i2
    
    call this%fft_r2c_sub( m, x )
    
    do concurrent ( i1 = 1:m )
      cx(i1,1)%re = x(i1,1,1)
      cx(i1,1)%im = zero
    end do
    
    do concurrent ( i2 = 2:this%np, i1 = 1:m )
      cx(i1,i2)%re = x(i1,1,i2)
      cx(i1,i2)%im = x(i1,2,i2)
    end do
    
  end subroutine fft_r2c_exec_sub
  
  pure module subroutine fft_r2c_sub(this, m, x)
    class(T_fft),   intent(in)    :: this
    integer,        intent(in)    :: m
    real(kind=dbl), intent(inout) :: x(m,2,0:this%n/2-1)
    integer                       :: iv, i, ii, isj, j, isj2
    real(kind=dbl)                :: scal, temp, addre, addim, subim, subre, tempre, tempim
    real(kind=dbl), allocatable   :: y(:,:)
    
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
      
    scal = 1._dbl / this%n
    
    do concurrent ( iv = 1:m )
       temp      =        x(iv,1,0) * scal
       x(iv,1,0) = temp + x(iv,2,0) * scal
       x(iv,2,0) = temp - x(iv,2,0) * scal
    end do
    
    do concurrent ( i = 1:(this%n-2)/4 , iv = 1:m )
      addre = ( x(iv,1,this%n/2-i) + x(iv,1,i) ) * scal / 2
      subre = ( x(iv,1,this%n/2-i) - x(iv,1,i) ) * scal / 2
      addim = ( x(iv,2,this%n/2-i) + x(iv,2,i) ) * scal / 2
      subim = ( x(iv,2,this%n/2-i) - x(iv,2,i) ) * scal / 2
      
      tempre = addre - subre * this%t(this%n+2*i) + addim * this%t(this%n+2*i-1)
      tempim = subim - addim * this%t(this%n+2*i) - subre * this%t(this%n+2*i-1)
      
      x(iv,1,         i) = tempre
      x(iv,2,         i) = tempim
      x(iv,1,this%n/2-i) = +2 * addre - tempre
      x(iv,2,this%n/2-i) = -2 * subim + tempim
    end do
    
    if ( mod(this%n,4) == 0) then
       do concurrent ( iv = 1:m )
         x(iv,1,this%n/4) = +x(iv,1,this%n/4) * scal
         x(iv,2,this%n/4) = -x(iv,2,this%n/4) * scal
       end do
    end if
    
  end subroutine fft_r2c_sub
  
end submodule fft_r2c