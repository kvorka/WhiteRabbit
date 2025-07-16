submodule (fourier_transform) r2c
  implicit none; contains
  
  module procedure fft_r2c_sub
    integer        :: i, iv, i1
    real(kind=dbl) :: addre, addim, subre, subim, t1, t2
    
    call this%fxztal(   m, x )
    call this%fxzshf( 2*m, x )
    
    do iv = 1, m
      !$omp simd
      do i1 = 1, 16
        addre        =                x(i1,iv,1,0)
        x(i1,iv,1,0) = x(i1,iv,1,0) + x(i1,iv,2,0)
        x(i1,iv,2,0) = addre        - x(i1,iv,2,0)
      end do
    end do
    
    do i = 1, (this%n-2)/4
      t1 = this%t(this%n+2*i-1)
      t2 = this%t(this%n+2*i  )
      
      do iv = 1, m
        !$omp simd
        do i1 = 1, 16
          addre = -x(i1,iv,1,i) - x(i1,iv,1,this%n/2-i)
          subre = +x(i1,iv,1,i) - x(i1,iv,1,this%n/2-i)
          addim = +x(i1,iv,2,i) + x(i1,iv,2,this%n/2-i)
          subim = -x(i1,iv,2,i) + x(i1,iv,2,this%n/2-i)
          
          x(i1,iv,1,i) = -( addre - subre * t2 - addim * t1 ) / 2
          x(i1,iv,2,i) = +( subim - addim * t2 + subre * t1 ) / 2
          
          x(i1,iv,1,this%n/2-i) = -x(i1,iv,1,i) - addre
          x(i1,iv,2,this%n/2-i) = +x(i1,iv,2,i) - subim
        end do
      end do
    end do
    
    if ( mod(this%n,4) == 0) then
      do iv = 1, m
        !$omp simd
        do i1 = 1, 16
          x(i1,iv,1,this%n/4) = +x(i1,iv,1,this%n/4)
          x(i1,iv,2,this%n/4) = -x(i1,iv,2,this%n/4)
        end do
      end do
    end if
    
  end procedure fft_r2c_sub
  
end submodule r2c