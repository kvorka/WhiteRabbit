submodule (fourier_transform) fxrsc
  implicit none; contains
  
  module pure subroutine fxzrsc(this, m, x, sgn, fac)
    class(T_fft),   intent(in)    :: this
    integer,        intent(in)    :: m, sgn
    real(kind=dbl), intent(in)    :: fac
    real(kind=dbl), intent(inout) :: x(m,2,0:this%n/2-1)
    integer                       :: i, iv, in
    real(kind=dbl)                :: tempre, tempim, addre, addim, subre, subim, fac2, t1, t2
    
    do concurrent ( iv = 1:m )
      tempre    =          x(iv,1,0)
      x(iv,1,0) = tempre + x(iv,2,0)
      x(iv,2,0) = tempre - x(iv,2,0)
    end do
    
    do i = 1, (this%n-2)/4
      in = this%n/2-i
      t1 = sgn * this%t(this%n+2*i-1)
      t2 = sgn * this%t(this%n+2*i  )
      
      do concurrent ( iv = 1:m )
        addre = ( x(iv,1,i) + x(iv,1,in) ) * fac
        subre = ( x(iv,1,i) - x(iv,1,in) ) * fac
        addim = ( x(iv,2,i) + x(iv,2,in) ) * fac
        subim = ( x(iv,2,i) - x(iv,2,in) ) * fac
        
        tempre = ( addre - subre * t2 - addim * t1 )
        tempim = ( subim - addim * t2 + subre * t1 )
        
        x(iv,1,i ) = +( tempre             )
        x(iv,2,i ) = +( tempim             ) * sgn
        x(iv,1,in) = -( tempre - 2 * addre )
        x(iv,2,in) = +( tempim - 2 * subim ) * sgn
      end do
    end do
    
    if ( mod(this%n,4) == 0) then
      fac2 = 2 * fac
      
      do concurrent ( iv = 1:m )
        x(iv,1,this%n/4) = +x(iv,1,this%n/4) * fac2
        x(iv,2,this%n/4) = -x(iv,2,this%n/4) * fac2
      end do
    end if
    
  end subroutine fxzrsc
  
end submodule fxrsc