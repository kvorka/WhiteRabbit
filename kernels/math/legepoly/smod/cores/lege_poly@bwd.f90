submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer :: m, j, ma, nb4
    
    ma  = 0
    nb4 = 4 * nb
    
    do m = 0, this%jmax
      call zero_rarray_sub( 4*nb*16, swork )
      
      !j = m
        ma = ma+1
        
        call mmset_sub( ma, this%fmj(1,ma), cosx, sinx, pmm, pmj2, pmj1, pmj )
        call bwd_sum_sub( nb4, pmj, cc(1,ma), swork )
      
      do j = 1, (this%jmax-m)/2
        ma = ma+1
        
        call mjrec_sub( this%fmj(1,ma), cosx2, pmj2, pmj1, pmj )
        call bwd_sum_sub( nb4, pmj, cc(1,ma), swork )
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        ma = ma+1
        
        call mjrec_sub( this%fmj(1,ma), cosx2, pmj2, pmj1, pmj )
        call bwd_sum_sub( nb4, pmj, cc(1,ma), swork )
      end if
      
      call bwd_shuffle_sub( nb, cosx, swork, sumN(1,m), sumS(1,m) )
    end do
    
  end procedure bwd_legesum_sub
  
end submodule bwd