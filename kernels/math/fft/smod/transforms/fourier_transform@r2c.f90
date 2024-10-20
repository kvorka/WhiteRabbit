submodule (fourier_transform) r2c
  implicit none; contains
  
  module pure subroutine fft_r2c_sub(this, m, x)
    class(T_fft),   intent(in)    :: this
    integer,        intent(in)    :: m
    real(kind=dbl), intent(inout) :: x(m,2,0:this%n/2-1)
    integer                       :: iv, i, isj, j, isj2
    real(kind=dbl)                :: scal, addre, addim, subim, subre, tempre, tempim, t1, t2
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
    
    scal = one / this%n
    
    do concurrent ( iv = 1:m )
       tempre    =          x(iv,1,0) * scal
       x(iv,1,0) = tempre + x(iv,2,0) * scal
       x(iv,2,0) = tempre - x(iv,2,0) * scal
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
  
end submodule r2c