submodule (fourier_transform) c2r
  implicit none ; contains
  
  module procedure fft_c2r_sub
    integer        :: i1, i2, i3
    real(kind=dbl) :: addre, addim, subre, subim, t1, t2
    
    do i2 = 1, m
      !$omp simd
      do i1 = 1, 16
        addre     =                   x(i1,i2,1,0)
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
          addre = x(i1,i2,1,i3) + x(i1,i2,1,this%n/2-i3)
          subre = x(i1,i2,1,i3) - x(i1,i2,1,this%n/2-i3)
          addim = x(i1,i2,2,i3) + x(i1,i2,2,this%n/2-i3)
          subim = x(i1,i2,2,i3) - x(i1,i2,2,this%n/2-i3)
          
          x(i1,i2,1,i3) = addre - subre * t2 - addim * t1
          x(i1,i2,2,i3) = subim - addim * t2 + subre * t1
          
          x(i1,i2,1,this%n/2-i3) = -x(i1,i2,1,i3) + 2 * addre
          x(i1,i2,2,this%n/2-i3) = +x(i1,i2,2,i3) - 2 * subim
        end do
      end do
    end do
    
    if ( mod(this%n,4) == 0) then
      do i2 = 1, m
        !$omp simd
        do i1 = 1, 16
          x(i1,i2,1,this%n/4) = +x(i1,i2,1,this%n/4) * 2
          x(i1,i2,2,this%n/4) = -x(i1,i2,2,this%n/4) * 2
        end do
      end do
    end if
    
    call this%fxztal(   m, x )
    call this%fxzshf( 2*m, x )
    
  end procedure fft_c2r_sub
  
end submodule c2r