submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_sub
    integer        :: i1, i2
    real(kind=dbl) :: rcr, icr
    
    do concurrent ( i2 = 1:n )
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
  
end submodule fwd_sum
