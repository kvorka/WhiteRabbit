submodule (lateral_grid) transform
  implicit none; contains
  
  module procedure transform_sub
    integer                     :: itheta
    type(c_ptr)                 :: c_work
    real(kind=dbl), pointer     :: work(:)
    real(kind=dbl), pointer     :: pmm(:), pmj2(:), pmj1(:), pmj(:)
    real(kind=dbl), pointer     :: cosx(:), sinx(:), cosx2(:), wght(:)
    real(kind=dbl), pointer     :: sumN(:), sumS(:), swork(:), grid(:)
    real(kind=dbl), allocatable :: rcr(:), rcc(:)
    
    !Prepare input and output arrays
    call this%lgp%alloc_rscal_sub( nb, rcc )
    call this%lgp%alloc_rscal_sub( nf, rcr )
    
    call this%lgp%index_bwd_sub( nb, cc, rcc )
    
    !Allocating memory
    call alloc_aligned1d_sub( 2*(2*(nb+1)+nb+nb*this%fourtrans%n)*step, c_work, work )
      
      pmm   => work(                                         1 :                                       step )
      pmj   => work(                                    step+1 :   2*                                  step )
      pmj1  => work(   2*                               step+1 :   3*                                  step )
      pmj2  => work(   3*                               step+1 :   4*                                  step )
      swork => work(   4*                               step+1 :   4*(nb+1)*                           step )
      sumN  => work(   4*(nb+1)*                        step+1 : ( 4*(nb+1)+     nb*this%fourtrans%n )*step )
      sumS  => work( ( 4*(nb+1)+  nb*this%fourtrans%n )*step+1 : ( 4*(nb+1)   +2*nb*this%fourtrans%n )*step )
      grid  => work( ( 4*(nb+1)+2*nb*this%fourtrans%n )*step+1 : ( 4*(nb+1)+nb+2*nb*this%fourtrans%n )*step )
    
    !Cycle over latitudes :: calculating step at once
    do itheta = 1, (this%lgp%nLege/step)*step, step
      cosx  => this%lgp%rw(itheta:itheta+step-1,1)
      sinx  => this%lgp%rw(itheta:itheta+step-1,2)
      cosx2 => this%lgp%rw(itheta:itheta+step-1,3)
      wght  => this%lgp%rw(itheta:itheta+step-1,4)
      
      call zero_rarray_sub( nb*step*this%fourtrans%n, sumN )
      call zero_rarray_sub( nb*step*this%fourtrans%n, sumS )
      
      call this%lgp%bwd_legesum_sub( nb, rcc, sumN, sumS, cosx, sinx, cosx2, pmm, pmj2, pmj1, pmj, swork )
      
      call this%fourtrans%fft_c2r_sub( step*nb, sumN )
      call this%fourtrans%fft_c2r_sub( step*nb, sumS )
      
      call grid_sub( this%fourtrans%n, sumS, grid )
      call grid_sub( this%fourtrans%n, sumN, grid )
      
      call this%fourtrans%fft_r2c_sub( step*nf, sumN )
      call this%fourtrans%fft_r2c_sub( step*nf, sumS )
      
      call this%lgp%fwd_legesum_sub( nf, sumN, sumS, rcr, cosx, sinx, cosx2, wght, pmm, pmj2, pmj1, pmj, swork )
    end do
    
    call this%lgp%index_fwd_sub( nf, cr, rcr )
    
    !Cleaning
    call free_aligned1d_sub( c_work, work )
    
    deallocate( rcc, rcr )
    
  end procedure transform_sub
  
end submodule transform