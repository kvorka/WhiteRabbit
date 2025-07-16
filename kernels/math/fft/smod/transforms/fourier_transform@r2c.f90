submodule (fourier_transform) r2c
  implicit none; contains
  
  module procedure fft_r2c_sub
    integer        :: i3, i2, i1
    real(kind=dbl) :: addre, addim, subre, subim, t1, t2
    
    call this%fxztal(   m, x )
    call this%fxzshf( 2*m, x )
    
    do i2 = 1, m
      !$omp simd
      do i1 = 1, 16
        addre        =                x(i1,i2,1,0)
        x(i1,i2,1,0) = x(i1,i2,1,0) + x(i1,i2,2,0)
        x(i1,i2,2,0) = addre        - x(i1,i2,2,0)
      end do
    end do
    
    do i3 = 1, (this%n-2)/4
      t1 = this%t(this%n+2*i3-1)
      t2 = this%t(this%n+2*i3  )
      
      do i2 = 1, m
        !$omp simd
        do i1 = 1, 16
          addre = -x(i1,i2,1,i3) - x(i1,i2,1,this%n/2-i3)
          subre = +x(i1,i2,1,i3) - x(i1,i2,1,this%n/2-i3)
          addim = +x(i1,i2,2,i3) + x(i1,i2,2,this%n/2-i3)
          subim = -x(i1,i2,2,i3) + x(i1,i2,2,this%n/2-i3)
          
          x(i1,i2,1,i3) = -( addre - subre * t2 - addim * t1 ) / 2
          x(i1,i2,2,i3) = +( subim - addim * t2 + subre * t1 ) / 2
          
          x(i1,i2,1,this%n/2-i3) = -x(i1,i2,1,i3) - addre
          x(i1,i2,2,this%n/2-i3) = +x(i1,i2,2,i3) - subim
        end do
      end do
    end do
    
    if ( mod(this%n,4) == 0) then
      do i2 = 1, m
        !$omp simd
        do i1 = 1, 16
          x(i1,i2,1,this%n/4) = +x(i1,i2,1,this%n/4)
          x(i1,i2,2,this%n/4) = -x(i1,i2,2,this%n/4)
        end do
      end do
    end if
    
  end procedure fft_r2c_sub
  
end submodule r2c