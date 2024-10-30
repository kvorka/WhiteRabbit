submodule (fourier_transform) c2r
  implicit none ; contains
  
  module procedure fft_c2r_sub
    
    call this%fxzrsc( m, x, 1, one )
    call this%fxztal( m, x )
    call this%fxzshf( m, x )
    
  end procedure fft_c2r_sub
  
end submodule c2r