submodule (lege_poly) fwd
  implicit none; contains

  module procedure fwd_legesum_sub
    integer                             :: j, m, ma, nf4
    real(kind=dbl), pointer, contiguous :: p2pmm(:), p2pmj1(:), p2pmj(:), p2tmp(:)
    
    ma  = 0
    nf4 = 4 * nf
    
    p2pmm  => pwork(:,1)
    p2pmj1 => pwork(:,2)
    p2pmj  => pwork(:,3)
    
    do m = 0, this%jmax
      call fwd_shuffle_sub( nf, weight, cosx, sumN(1,m), sumS(1,m), swork )
      
      !j = m
        ma = ma+1
        
        call mmset_sub( ma, this%fmj(2,ma), cosx, sinx, p2pmm, p2pmj1, p2pmj )
        call fwd_sum_sub( nf4, p2pmj, swork, cr(1,ma) )
      
      do j = 1, (this%jmax-m)/2
        ma = ma+1
        
        p2tmp  => p2pmj1
        p2pmj1 => p2pmj
        p2pmj  => p2tmp
        
        call mjrec_sub( this%fmj(1,ma), cosx2, p2pmj1, p2pmj )
        call fwd_sum_sub( nf4, p2pmj, swork, cr(1,ma) )
      end do
      
      if ( mod(this%jmax-m,2) /= 0 ) then
        ma = ma+1
        
        p2tmp  => p2pmj1
        p2pmj1 => p2pmj
        p2pmj  => p2tmp
        
        call mjrec_sub( this%fmj(1,ma), cosx2, p2pmj1, p2pmj )
        call fwd_sum_sub( nf4, p2pmj, swork, cr(1,ma) )
      end if
    end do
    
  end procedure fwd_legesum_sub
  
end submodule fwd