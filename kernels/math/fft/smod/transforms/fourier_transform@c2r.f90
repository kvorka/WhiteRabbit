submodule (fourier_transform) c2r
  implicit none ; contains
  
  module pure subroutine fft_c2r_sub(this, m, x)
    class(T_fft),   intent(in)    :: this
    integer,        intent(in)    :: m
    real(kind=dbl), intent(inout) :: x(*)
    
    call this%fxzrsc( m, x, 1, one )
    call this%fxztal( m, x )
    call this%fxzshf( m, x )
    
  end subroutine fft_c2r_sub
  
end submodule c2r