submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: m, j, ma, nb4
    real(kind=dbl), pointer, contiguous :: p2pmm(:), p2pmj1(:), p2pmj(:), p2tmp(:)
    
    ma  = 0
    nb4 = 4 * nb
    
    p2pmm  => pwork(:,1)
    p2pmj1 => pwork(:,2)
    p2pmj  => pwork(:,3)
    
    do m = 0, this%jmax
      call zero_rarray_sub( 4*nb*step, swork )
      
      !j = m
        ma = ma+1
        
        call mmset_sub( ma, this%fmj(2,ma), cosx, sinx, p2pmm, p2pmj1, p2pmj )
        call bwd_sum_sub( nb4, p2pmj, cc(1,ma), swork )
      
      do j = 1, (this%jmax-m)/2
        ma = ma+1
        
        p2tmp  => p2pmj1
        p2pmj1 => p2pmj
        p2pmj  => p2tmp
        
        call mjrec_sub( this%fmj(1,ma), cosx2, p2pmj1, p2pmj )
        call bwd_sum_sub( nb4, p2pmj, cc(1,ma), swork )
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        ma = ma+1
        
        p2tmp  => p2pmj1
        p2pmj1 => p2pmj
        p2pmj  => p2tmp
        
        call mjrec_sub( this%fmj(1,ma), cosx2, p2pmj1, p2pmj )
        call bwd_sum_sub( nb4, p2pmj, cc(1,ma), swork )
      end if
      
      call bwd_shuffle_sub( nb, cosx, swork, sumN(1,m), sumS(1,m) )
    end do
    
  end procedure bwd_legesum_sub
  
end submodule bwd