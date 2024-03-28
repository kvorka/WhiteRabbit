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
    real(kind=dbl),    intent(in)  :: cx(2,m,*)
    real(kind=dbl),    intent(out) :: x(m,2,*)
    integer                        :: i1, i2, i3
    
    do concurrent ( i3 = 1:this%np, i1 = 1:m )
      x(i1,1,i3) = cx(1,i1,i3)
      x(i1,2,i3) = cx(2,i1,i3)
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
    real(kind=dbl),    intent(out)   :: cx(2,m,*)
    integer                          :: i1, i2
    
    call this%fft_r2c_sub( m, x )
    
    do concurrent ( i1 = 1:m )
      cx(1,i1,1) = x(i1,1,1)
      cx(2,i1,1) = zero
    end do
    
    do concurrent ( i2 = 2:this%np, i1 = 1:m )
      cx(1,i1,i2) = x(i1,1,i2)
      cx(2,i1,i2) = x(i1,2,i2)
    end do
    
  end subroutine fft_r2c_exec_sub
  
  module pure subroutine fft_deallocate_sub(this)
    class(T_fft), intent(inout) :: this
    
    if ( allocated( this%it ) ) deallocate( this%it )
    if ( allocated( this%t  ) ) deallocate( this%t  )
    
  end subroutine fft_deallocate_sub
  
end submodule fft_subs