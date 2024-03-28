submodule (Fourier_transform) Init_FourierTransform
  implicit none; contains
  
  module pure subroutine fft_init_sub(this, n)
    class(T_fft), intent(inout) :: this
    integer,      intent(in)    :: n
    integer                     :: i
    
    this%n = n
      allocate( this%it(n/2)  ) ; this%it = 0
      allocate( this%t(3*n/2) ) ; this%t = zero
    
    call fxzini(n/2, this%it, this%t)
    
    do i = 1, (n-2) / 4
      this%t(n+2*i-1) = cos(2 * pi * i / n)
      this%t(n+2*i  ) = sin(2 * pi * i / n)        
    end do
    
  end subroutine fft_init_sub
  
  module pure subroutine fft_deallocate_sub(this)
    class(T_fft), intent(inout) :: this
    
    if ( allocated( this%it ) ) deallocate( this%it )
    if ( allocated( this%t  ) ) deallocate( this%t  )
    
  end subroutine fft_deallocate_sub
  
end submodule Init_FourierTransform