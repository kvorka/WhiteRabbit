submodule (lateral_grid) transform
  implicit none; contains
  
  module procedure transform_sub
    integer                             :: itheta
    type(c_ptr)                         :: c_work
    real(kind=dbl), pointer, contiguous :: work(:)
    real(kind=dbl), pointer, contiguous :: pmm(:), pmj2(:), pmj1(:), pmj(:)
    real(kind=dbl), pointer, contiguous :: cosx(:), sinx(:), cosx2(:), wght(:)
    real(kind=dbl), pointer, contiguous :: sumN(:), sumS(:), swork(:), grid(:)
    real(kind=dbl), allocatable         :: rcr(:), rcc(:)
    
    !Prepare input and output arrays
    call this%lgp%alloc_rscal_sub( nb, rcc )
    call this%lgp%alloc_rscal_sub( nf, rcr )
    
    call this%lgp%index_bwd_sub( nb, cc, rcc )
    
    !Allocating memory
    call alloc_aligned1d_sub( 2*(2*(nb+1)+nb+nb*this%fourtrans%n)*16, c_work, work )
      
      pmm   => work(                                       1 :                                       16 )
      pmj   => work(                                    16+1 :   2*                                  16 )
      pmj1  => work(   2*                               16+1 :   3*                                  16 )
      pmj2  => work(   3*                               16+1 :   4*                                  16 )
      swork => work(   4*                               16+1 :   4*(nb+1)*                           16 )
      sumN  => work(   4*(nb+1)*                        16+1 : ( 4*(nb+1)+     nb*this%fourtrans%n )*16 )
      sumS  => work( ( 4*(nb+1)+  nb*this%fourtrans%n )*16+1 : ( 4*(nb+1)   +2*nb*this%fourtrans%n )*16 )
      grid  => work( ( 4*(nb+1)+2*nb*this%fourtrans%n )*16+1 : ( 4*(nb+1)+nb+2*nb*this%fourtrans%n )*16 )
    
    !Cycle over latitudes :: calculating 16 at once
    do itheta = 1, (this%lgp%nLege/16)*16, 16
      cosx  => this%lgp%rw(itheta:itheta+15,1)
      sinx  => this%lgp%rw(itheta:itheta+15,2)
      cosx2 => this%lgp%rw(itheta:itheta+15,3)
      wght  => this%lgp%rw(itheta:itheta+15,4)
      
      call zero_rarray_sub( nb*16*this%fourtrans%n, sumN )
      call zero_rarray_sub( nb*16*this%fourtrans%n, sumS )
      
      call this%lgp%bwd_legesum_sub( nb, rcc, sumN, sumS, cosx, sinx, cosx2, pmm, pmj2, pmj1, pmj, swork )
      
      call this%fourtrans%fft_c2r_sub( nb, sumN )
      call this%fourtrans%fft_c2r_sub( nb, sumS )
      
      call grid_sub( this%fourtrans%n, sumS, grid )
      call grid_sub( this%fourtrans%n, sumN, grid )
      
      call this%fourtrans%fft_r2c_sub( nf, sumN )
      call this%fourtrans%fft_r2c_sub( nf, sumS )
      
      call this%lgp%fwd_legesum_sub( nf, sumN, sumS, rcr, cosx, sinx, cosx2, wght, pmm, pmj2, pmj1, pmj, swork )
    end do
    
    call this%lgp%index_fwd_sub( nf, cr, rcr )
    
    !Cleaning
    call free_aligned1d_sub( c_work, work )
    
    deallocate( rcc, rcr )
    
  end procedure transform_sub
  
end submodule transform