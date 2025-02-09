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
    c_swork = malloc( alig, 4*nb*size_step                  )
    c_pmm   = malloc( alig,      size_step                  )
    c_pmj   = malloc( alig,      size_step                  )
    c_pmj1  = malloc( alig,      size_step                  )
    c_pmj2  = malloc( alig,      size_step                  )
    c_cosx  = malloc( alig,      size_step                  )
    c_cosx2 = malloc( alig,      size_step                  )
    c_sinx  = malloc( alig,      size_step                  )
    c_wght  = malloc( alig,      size_step                  )
    c_sumN  = malloc( alig,   nb*size_step*this%fourtrans%n )
    c_sumS  = malloc( alig,   nb*size_step*this%fourtrans%n )
    
    call c_f_pointer( c_swork, swork, [ 4*nb*step                  ] )
    call c_f_pointer( c_pmm,   pmm,   [      step                  ] )
    call c_f_pointer( c_pmj,   pmj,   [      step                  ] )
    call c_f_pointer( c_pmj1,  pmj1,  [      step                  ] )
    call c_f_pointer( c_pmj2,  pmj2,  [      step                  ] )
    call c_f_pointer( c_cosx,  cosx,  [      step                  ] )
    call c_f_pointer( c_cosx2, cosx2, [      step                  ] )
    call c_f_pointer( c_sinx,  sinx,  [      step                  ] )
    call c_f_pointer( c_wght,  wght,  [      step                  ] )
    call c_f_pointer( c_sumN,  sumN,  [   nb*step*this%fourtrans%n ] )
    call c_f_pointer( c_sumS,  sumS,  [   nb*step*this%fourtrans%n ] )
    
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
    nullify( swork )
    nullify( pmm   )
    nullify( pmj   )
    nullify( pmj1  )
    nullify( pmj2  )
    nullify( cosx  )
    nullify( sinx  )
    nullify( cosx2 )
    nullify( wght  )
    nullify( sumN  )
    nullify( sumS  )
    
    call free( c_swork )
    call free( c_pmm   )
    call free( c_pmj2  )
    call free( c_pmj1  )
    call free( c_pmj   )
    call free( c_cosx  )
    call free( c_sumN  )
    call free( c_sumS  )
    call free( c_sinx  )
    call free( c_cosx2 )
    call free( c_wght  )
    
    deallocate( rcc, rcr )
    
  end procedure transform_sub
  
end submodule transform