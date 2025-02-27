submodule (lateral_grid) transform
  implicit none; contains
  
  module procedure transform_sub
    integer                             :: itheta
    type(c_ptr)                         :: c_work, c_rcc, c_rcr
    real(kind=dbl), contiguous, pointer :: work(:)
    real(kind=dbl), contiguous, pointer :: pwork(:), grid(:)
    real(kind=dbl), contiguous, pointer :: cosx(:), sinx(:), cosx2(:), wght(:)
    real(kind=dbl), contiguous, pointer :: sumN(:), sumS(:), swork(:)
    real(kind=dbl), contiguous, pointer :: rcc(:), rcr(:)
    
    !Prepare input and output arrays
    call this%lgp%alloc_rscal_sub( nb, c_rcc, rcc )
    call this%lgp%alloc_rscal_sub( nf, c_rcr, rcr )
    
    call this%lgp%index_bwd_sub( nb, cc, rcc )
    
    !Allocating memory
    call alloc_aligned1d_sub( (2*nb*this%fourtrans%n+5*nb+3)*step, c_work, work )
    
    pwork => work(                                           1 : (                              3 ) * step )
    swork => work( (                            3 ) * step + 1 : (                         4*nb+3 ) * step )
    grid  => work( (                       4*nb+3 ) * step + 1 : (                         5*nb+3 ) * step )
    sumN  => work( (                       5*nb+3 ) * step + 1 : (   nb*this%fourtrans%n + 5*nb+3 ) * step )
    sumS  => work( ( nb*this%fourtrans%n + 5*nb+3 ) * step + 1 : ( 2*nb*this%fourtrans%n + 5*nb+3 ) * step )
      
    !Cycle over latitudes :: calculating step at once
    do itheta = 1, (this%lgp%nLege/step)*step, step
      cosx  => this%lgp%rw(itheta:itheta+step-1,1)
      sinx  => this%lgp%rw(itheta:itheta+step-1,2)
      cosx2 => this%lgp%rw(itheta:itheta+step-1,3)
      wght  => this%lgp%rw(itheta:itheta+step-1,4)
      
      call zero_rarray_sub( nb*step*this%fourtrans%n, sumN )
      call zero_rarray_sub( nb*step*this%fourtrans%n, sumS )
      
      call this%lgp%bwd_legesum_sub( nb, rcc, sumN, sumS, cosx, sinx, cosx2, pwork, swork )
      
      call this%fourtrans%fft_c2r_sub( step*nb, sumN )
      call this%fourtrans%fft_c2r_sub( step*nb, sumS )
      
      call grid_sub( this%fourtrans%n, sumS, grid )
      call grid_sub( this%fourtrans%n, sumN, grid )
      
      call this%fourtrans%fft_r2c_sub( step*nf, sumN )
      call this%fourtrans%fft_r2c_sub( step*nf, sumS )
      
      call this%lgp%fwd_legesum_sub( nf, sumN, sumS, rcr, cosx, sinx, cosx2, wght, pwork, swork )
    end do
    
    call this%lgp%index_fwd_sub( nf, cr, rcr )
    
    !Cleaning
    call free_aligned1d_sub( c_work, work )
    call free_aligned1d_sub( c_rcc, rcc )
    call free_aligned1d_sub( c_rcr, rcr )
    
  end procedure transform_sub
  
end submodule transform