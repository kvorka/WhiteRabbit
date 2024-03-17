submodule (Fourier_transform) fft_subs
  implicit none; contains
  
  module pure subroutine fft_init_sub(this, n)
    class(T_fft), intent(inout) :: this
    integer,      intent(in)    :: n
    integer                     :: i
    
    this%np = n/3
    this%n  = n
      allocate( this%it(n/2)  ) ; this%it = 0
      allocate( this%t(3*n/2) ) ; this%t = 0._dbl
    
    call fxzini(n/2, this%it, this%t)
    
    do i = 1, (n-2) / 4
      this%t(n+2*i-1) = cos(2 * pi * i / n)
      this%t(n+2*i  ) = sin(2 * pi * i / n)        
    end do
    
  end subroutine fft_init_sub
  
  module pure subroutine fft_c2r_exec_sub(this, m, cx, x)
    class(T_fft),      intent(in)  :: this
    integer,           intent(in)  :: m
    complex(kind=dbl), intent(in)  :: cx(m,*)
    real(kind=dbl),    intent(out) :: x(m,2,*)
    integer                        :: i1, i2, i3
    
    do concurrent ( i3 = 1:this%np, i1 = 1:m )
      x(i1,1,i3) = cx(i1,i3)%re
      x(i1,2,i3) = cx(i1,i3)%im
    end do
    
    do concurrent ( i3 = this%np+1:this%n/2, i2 = 1:2, i1 = 1:m )
      x(i1,i2,i3) = zero
    end do
    
    call this%fft_c2r_sub( m,  x )
    
  end subroutine fft_c2r_exec_sub
  
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
  
  module pure subroutine fft_deallocate_sub(this)
    class(T_fft), intent(inout) :: this
    
    if ( allocated( this%it ) ) deallocate( this%it )
    if ( allocated( this%t  ) ) deallocate( this%t  )
    
  end subroutine fft_deallocate_sub
  
end submodule fft_subs