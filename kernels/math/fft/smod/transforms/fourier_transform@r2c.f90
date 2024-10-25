submodule (fourier_transform) r2c
  implicit none; contains
  
  module pure subroutine fft_r2c_sub(this, m, x)
    class(T_fft),   intent(in)    :: this
    integer,        intent(in)    :: m
    real(kind=dbl), intent(inout) :: x(*)
    
    call this%fxztal( m, x )
    call this%fxzshf( m, x )
    call this%fxzrsc( m, x, -1, 0.5_dbl )
    
  end subroutine fft_r2c_sub
  
end submodule r2c