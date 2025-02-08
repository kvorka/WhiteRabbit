submodule (lateral_grid) transform
  implicit none; contains
  
  module procedure transform_sub
    integer                        :: itheta
    real(kind=dbl),    pointer     :: swork(:)
    real(kind=dbl),    pointer     :: pmm(:), pmj2(:), pmj1(:), pmj(:)
    real(kind=dbl),    pointer     :: cosx(:), sinx(:), cosx2(:), wght(:)
    real(kind=dbl),    pointer     :: sumN(:), sumS(:)
    complex(kind=dbl), allocatable :: rcr(:), rcc(:)
    type(c_ptr)                    :: c_swork, c_pmm, c_pmj2, c_pmj1, c_pmj, c_cosx, c_sinx, c_cosx2, c_wght, c_sumN, c_sumS
    
    !Prepare input and output arrays
    call this%lgp%alloc_cscal_sub( nb, rcc )
    call this%lgp%alloc_cscal_sub( nf, rcr )
    
    call this%lgp%index_bwd_sub( nb, cc, rcc )
    
    !Allocating memory
    call calloc_sub( 4*nb*step,                  step, c_swork ); call c_f_pointer( c_swork, swork, [ 4*nb*step                  ] )
    call calloc_sub(      step,                  step, c_pmm   ); call c_f_pointer( c_pmm,   pmm,   [      step                  ] )
    call calloc_sub(      step,                  step, c_pmj   ); call c_f_pointer( c_pmj,   pmj,   [      step                  ] )
    call calloc_sub(      step,                  step, c_pmj1  ); call c_f_pointer( c_pmj1,  pmj1,  [      step                  ] )
    call calloc_sub(      step,                  step, c_pmj2  ); call c_f_pointer( c_pmj2,  pmj2,  [      step                  ] )
    call calloc_sub(      step,                  step, c_cosx  ); call c_f_pointer( c_cosx,  cosx,  [      step                  ] )
    call calloc_sub(      step,                  step, c_cosx2 ); call c_f_pointer( c_cosx2, cosx2, [      step                  ] )
    call calloc_sub(      step,                  step, c_sinx  ); call c_f_pointer( c_sinx,  sinx,  [      step                  ] )
    call calloc_sub(      step,                  step, c_wght  ); call c_f_pointer( c_wght,  wght,  [      step                  ] )
    call calloc_sub(   nb*step*this%fourtrans%n, step, c_sumN  ); call c_f_pointer( c_sumN,  sumN,  [   nb*step*this%fourtrans%n ] )
    call calloc_sub(   nb*step*this%fourtrans%n, step, c_sumS  ); call c_f_pointer( c_sumS,  sumS,  [   nb*step*this%fourtrans%n ] )
    
    !Cycle over latitudes :: calculating step at once
    do itheta = 1, (this%lgp%nLege/step)*step, step
      cosx  = this%lgp%rw(itheta:itheta+step-1,1)
      sinx  = this%lgp%rw(itheta:itheta+step-1,2)
      cosx2 = this%lgp%rw(itheta:itheta+step-1,3)
      wght  = this%lgp%rw(itheta:itheta+step-1,4)
      
      call czero_sub( nb*step*this%fourtrans%n, sumN )
      call czero_sub( nb*step*this%fourtrans%n, sumS )
      
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
    call cfree_sub( c_swork )
    call cfree_sub( c_pmm   )
    call cfree_sub( c_pmj2  )
    call cfree_sub( c_pmj1  )
    call cfree_sub( c_pmj   )
    call cfree_sub( c_cosx  )
    call cfree_sub( c_sumN  )
    call cfree_sub( c_sumS  )
    call cfree_sub( c_sinx  )
    call cfree_sub( c_cosx2 )
    call cfree_sub( c_wght  )
    
    deallocate( rcc, rcr )
    
  end procedure transform_sub
  
end submodule transform