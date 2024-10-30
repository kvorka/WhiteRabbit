submodule (fourier_transform) r2c
  implicit none; contains
  
  module procedure fft_r2c_sub
    
    call this%fxztal( m, x )
    call this%fxzshf( m, x )
    call this%fxzrsc( m, x, -1, 0.5_dbl )
    
  end procedure fft_r2c_sub
  
end submodule r2c