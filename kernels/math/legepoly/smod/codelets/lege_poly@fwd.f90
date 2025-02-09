submodule (lege_poly) fwd
  implicit none; contains
  
  module procedure fwd_sum_sub
    integer        :: i1, i2
    real(kind=dbl) :: rcr, icr
    
    do i2 = 1, n
      rcr = zero
      icr = zero
      
      !$omp simd
      do i1 = 1, step
        rcr = rcr + pmj(i1) * swork(i1,1,i2)
        icr = icr + pmj(i1) * swork(i1,2,i2)
      end do
      
      cr(i2) = cr(i2) + cmplx( rcr, icr, kind=dbl )
    end do
    
  end procedure fwd_sum_sub
  
  module procedure fwd_shuffle_sub
    integer :: i1, i2, i3, i4
    
    do i2 = 1, n
      !$omp simd
      do i1 = 1, step
        swork(i1,1,i2,1) = sumN(i1,i2,1) - sumS(i1,i2,1)
        swork(i1,2,i2,1) = sumN(i1,i2,2) - sumS(i1,i2,2)
        swork(i1,1,i2,2) = sumN(i1,i2,1) + sumS(i1,i2,1)
        swork(i1,2,i2,2) = sumN(i1,i2,2) + sumS(i1,i2,2)
      end do
    end do
    
    !$omp simd collapse (3)
    do i3 = 1, n
      do i2 = 1, 2
        do i1 = 1, step
          swork(i1,i2,i3,1) = swork(i1,i2,i3,1) * w(i1)
          swork(i1,i2,i3,2) = swork(i1,i2,i3,2) * w(i1) * cosx(i1)
        end do
      end do
    end do
    
  end procedure fwd_shuffle_sub
  
  module procedure fwd_legesum_sub
    integer :: j, m, ma, nf2
    
    ma  = 0
    nf2 = 2 * nf
    
    do m = 0, this%jmax
      call fwd_shuffle_sub( nf, weight, cosx, sumN(1,m), sumS(1,m), swork )
      
      !j = m
        ma = ma+1
        
        call mmset_sub( ma, this%fmj(1,ma), cosx, sinx, pmm, pmj2, pmj1, pmj )
        call fwd_sum_sub( nf2, pmj, swork, cr(1,ma) )
      
      do j = 1, (this%jmax-m)/2
        ma = ma+1
        
        call mjrec_sub( this%fmj(1,ma), cosx2, pmj2, pmj1, pmj )
        call fwd_sum_sub( nf2, pmj, swork, cr(1,ma) )
      end do
      
      if ( mod(this%jmax-m,2) /= 0 ) then
        ma = ma+1
        
        call mjrec_sub( this%fmj(1,ma), cosx2, pmj2, pmj1, pmj )
        call fwd_sum_sub( nf2, pmj, swork, cr(1,ma) )
      end if
    end do
    
  end procedure fwd_legesum_sub
  
end submodule fwd