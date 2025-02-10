submodule (lateral_grid) transform
  implicit none; contains
  
  module procedure transform_sub
    integer                        :: itheta
    real(kind=dbl),    pointer     :: swork(:), pmm(:), pmj2(:), pmj1(:), pmj(:)
    real(kind=dbl),    pointer     :: cosx(:), sinx(:), cosx2(:), wght(:)
    real(kind=dbl),    pointer     :: sumN(:), sumS(:)
    complex(kind=dbl), allocatable :: rcr(:), rcc(:)
    type(c_ptr)                    :: c_swork, c_pmm, c_pmj2, c_pmj1, c_pmj
    type(c_ptr)                    :: c_cosx, c_sinx, c_cosx2, c_wght
    type(c_ptr)                    :: c_sumN, c_sumS
    
    !Prepare input and output arrays
    call this%lgp%alloc_cscal_sub( nb, rcc )
    call this%lgp%alloc_cscal_sub( nf, rcr )
    
    call this%lgp%index_bwd_sub( nb, cc, rcc )
    
    !Allocating memory
    call alloc_aligned_sub( 4*nb*step,                  c_swork, swork )
    call alloc_aligned_sub(      step,                  c_pmm,   pmm   )
    call alloc_aligned_sub(      step,                  c_pmj,   pmj   )
    call alloc_aligned_sub(      step,                  c_pmj1,  pmj1  )
    call alloc_aligned_sub(      step,                  c_pmj2,  pmj2  )
    call alloc_aligned_sub(      step,                  c_cosx,  cosx  )
    call alloc_aligned_sub(      step,                  c_cosx2, cosx2 )
    call alloc_aligned_sub(      step,                  c_sinx,  sinx  )
    call alloc_aligned_sub(      step,                  c_wght,  wght  )
    call alloc_aligned_sub(   nb*step*this%fourtrans%n, c_sumN,  sumN  )
    call alloc_aligned_sub(   nb*step*this%fourtrans%n, c_sumS,  sumS  )
    
    !Cycle over latitudes :: calculating step at once
    do itheta = 1, (this%lgp%nLege/step)*step, step
      cosx  = this%lgp%rw(itheta:itheta+step-1,1)
      sinx  = this%lgp%rw(itheta:itheta+step-1,2)
      cosx2 = this%lgp%rw(itheta:itheta+step-1,3)
      wght  = this%lgp%rw(itheta:itheta+step-1,4)
      
      call zero_rarray_sub( nb*step*this%fourtrans%n, sumN )
      call zero_rarray_sub( nb*step*this%fourtrans%n, sumS )
      
      call this%lgp%bwd_legesum_sub( nb, rcc, sumN, sumS, cosx, sinx, cosx2, pmm, pmj2, pmj1, pmj, swork )
      
      call this%fourtrans%fft_c2r_sub( step*nb, sumN )
      call this%fourtrans%fft_c2r_sub( step*nb, sumS )
      
      call grid_sub( this%fourtrans%n, sumS )
      call grid_sub( this%fourtrans%n, sumN )
      
      call this%fourtrans%fft_r2c_sub( step*nf, sumN )
      call this%fourtrans%fft_r2c_sub( step*nf, sumS )
      
      call this%lgp%fwd_legesum_sub( nf, sumN, sumS, rcr, cosx, sinx, cosx2, wght, pmm, pmj2, pmj1, pmj, swork )
    end do
    
    call this%lgp%index_fwd_sub( nf, cr, rcr )
    
    !Cleaning
    call free_aligned_sub( c_swork, swork )
    call free_aligned_sub( c_pmm,   pmm   )
    call free_aligned_sub( c_pmj,   pmj   )
    call free_aligned_sub( c_pmj1,  pmj1  )
    call free_aligned_sub( c_pmj2,  pmj2  )
    call free_aligned_sub( c_cosx,  cosx  )
    call free_aligned_sub( c_sinx,  sinx  )
    call free_aligned_sub( c_cosx2, cosx2 )
    call free_aligned_sub( c_wght,  wght  )
    call free_aligned_sub( c_sumN,  sumN  )
    call free_aligned_sub( c_sumS,  sumS  )
    
    deallocate( rcc, rcr )
    
  end procedure transform_sub
  
end submodule transform