submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_sum_sub
    integer        :: i1, i2
    real(kind=dbl) :: rcc, icc
    
    do i2 = 1, n
      rcc = cc(i2)%re
      icc = cc(i2)%im
      
      !$omp simd
      do i1 = 1, step
        swork(i1,1,i2) = swork(i1,1,i2) + pmj(i1) * rcc
        swork(i1,2,i2) = swork(i1,2,i2) + pmj(i1) * icc
      end do
    end do
    
  end procedure bwd_sum_sub
  
  module procedure bwd_shuffle_sub
    integer :: i1, i2
    
    do i2 = 1, n
      !$omp simd
      do i1 = 1, step
        swork(i1,1,i2,2) = swork(i1,1,i2,2) * cosx(i1)
        swork(i1,2,i2,2) = swork(i1,2,i2,2) * cosx(i1)
      end do
      
      !$omp simd
      do i1 = 1, step
        sumN(i1,i2,1) = swork(i1,1,i2,2) + swork(i1,1,i2,1)
        sumN(i1,i2,2) = swork(i1,2,i2,2) + swork(i1,2,i2,1)
        sumS(i1,i2,1) = swork(i1,1,i2,2) - swork(i1,1,i2,1)
        sumS(i1,i2,2) = swork(i1,2,i2,2) - swork(i1,2,i2,1)
      end do
    end do
    
  end procedure bwd_shuffle_sub
  
  module procedure bwd_legesum_sub
    integer :: m, j, ma, nb2
    
    ma  = 0
    nb2 = 2 * nb
    
    do m = 0, this%jmax
      call zero_rarray_sub( 4*nb*step, swork )
      
      !j = m
        ma = ma+1
        
        call mmset_sub( ma, this%fmj(1,ma), cosx, sinx, pmm, pmj2, pmj1, pmj )
        call bwd_sum_sub( nb2, pmj, cc(1,ma), swork )
      
      do j = 1, (this%jmax-m)/2
        ma = ma+1
        
        call mjrec_sub( this%fmj(1,ma), cosx2, pmj2, pmj1, pmj )
        call bwd_sum_sub( nb2, pmj, cc(1,ma), swork )
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        ma = ma+1
        
        call mjrec_sub( this%fmj(1,ma), cosx2, pmj2, pmj1, pmj )
        call bwd_sum_sub( nb2, pmj, cc(1,ma), swork )
      end if
      
      call bwd_shuffle_sub( nb, cosx, swork, sumN(1,m), sumS(1,m) )
    end do
    
  end procedure bwd_legesum_sub

end submodule bwd