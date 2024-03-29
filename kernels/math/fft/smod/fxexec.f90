submodule (Fourier_transform) fxexec
  implicit none ; contains
  
  module pure subroutine fft_c2r_exec_sub(this, m, cx, x)
    class(T_fft),      intent(in)  :: this
    integer,           intent(in)  :: m
    complex(kind=dbl), intent(in)  :: cx(m,*)
    real(kind=dbl),    intent(out) :: x(m,2,*)
    integer                        :: i1, i2, i3
    
    do concurrent ( i1 = 1:m )
      x(i1,1,1) = cx(i1,1)%re
      x(i1,2,1) = zero
    end do
    
    do concurrent ( i3 = 2:this%np, i1 = 1:m )
      x(i1,1,i3) = cx(i1,i3)%re
      x(i1,2,i3) = cx(i1,i3)%im
    end do
    
    do concurrent ( i3 = this%np+1:this%n/2, i2 = 1:2, i1 = 1:m )
      x(i1,i2,i3) = zero
    end do
    
    call this%fft_c2r_sub( m,  x )
    
  end subroutine fft_c2r_exec_sub
  
  module pure subroutine fft_c2r_sub(this, m, x)
    class(T_fft),   intent(in)    :: this
    integer,        intent(in)    :: m
    real(kind=dbl), intent(inout) :: x(m,2,0:this%n/2-1)
    integer                       :: iv, i, isj, j, isj2
    real(kind=dbl)                :: temp, tempre1, tempim1, tempre2, tempim2, t1, t2
    real(kind=dbl), allocatable   :: y(:,:)
    
    do concurrent ( iv = 1:m )
      temp      = x(iv,1,0)
      x(iv,1,0) = x(iv,1,0) + x(iv,2,0)
      x(iv,2,0) = temp      - x(iv,2,0)
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
        y(1:m,1:2) = x(1:m,1:2,isj)
        
        do
          j    = j + 1
          isj2 = this%it(j)
          
          if ( isj2 < 0 ) then
            x(1:m,1:2,isj)      = x(1:m,1:2,isj2-imm)
            x(1:m,1:2,isj2-imm) = y(1:m,1:2)
            
            j = j + 1 ; exit
          else
            x(1:m,1:2,isj) = x(1:m,1:2,isj2)
            isj = isj2
          end if
        end do
      end if
    end do
    
    deallocate( y )
    
  end subroutine fft_c2r_sub
  
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
  
  module pure subroutine fft_r2c_sub(this, m, x)
    class(T_fft),   intent(in)    :: this
    integer,        intent(in)    :: m
    real(kind=dbl), intent(inout) :: x(m,2,0:this%n/2-1)
    integer                       :: iv, i, isj, j, isj2
    real(kind=dbl)                :: scal, temp, addre, addim, subim, subre, tempre, tempim, t1, t2
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
        y(1:m,1:2) = x(1:m,1:2,isj)
        
        do
          j    = j + 1
          isj2 = this%it(j)
          
          if ( isj2 < 0 ) then
            x(1:m,1:2,isj)      = x(1:m,1:2,isj2-imm)
            x(1:m,1:2,isj2-imm) = y(1:m,1:2)
            
            j = j + 1 ; exit
          else
            x(1:m,1:2,isj) = x(1:m,1:2,isj2)
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
    
    scal = scal / 2
    
    do i = 1, (this%n-2)/4
      t1 = this%t(this%n+2*i-1)
      t2 = this%t(this%n+2*i  )
      
      do concurrent ( iv = 1:m )
        addre = ( x(iv,1,this%n/2-i) + x(iv,1,i) ) * scal
        subre = ( x(iv,1,this%n/2-i) - x(iv,1,i) ) * scal
        addim = ( x(iv,2,this%n/2-i) + x(iv,2,i) ) * scal
        subim = ( x(iv,2,this%n/2-i) - x(iv,2,i) ) * scal
        
        tempre = addre - subre * t2 + addim * t1
        tempim = subim - addim * t2 - subre * t1
        
        x(iv,1,         i) = tempre
        x(iv,2,         i) = tempim
        x(iv,1,this%n/2-i) = +2 * addre - tempre
        x(iv,2,this%n/2-i) = -2 * subim + tempim
      end do
    end do
    
    scal = scal * 2
    
    if ( mod(this%n,4) == 0) then
       do concurrent ( iv = 1:m )
         x(iv,1,this%n/4) = +x(iv,1,this%n/4) * scal
         x(iv,2,this%n/4) = -x(iv,2,this%n/4) * scal
       end do
    end if
    
  end subroutine fft_r2c_sub
  
end submodule fxexec