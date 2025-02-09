submodule (lateral_grid) sgs
  implicit none; contains
  
  module procedure allocate_grid_sub
    
    if ( allocated(grid) ) deallocate( grid )
    
    allocate( grid(this%lgp%nLege,this%fourtrans%n,2) )
      call zero_rarray_sub( 2 * this%lgp%nLege * this%fourtrans%n, grid )
    
  end procedure allocate_grid_sub
  
  module procedure space_to_grid_sub
    integer                        :: itheta, i1, i2
    real(kind=dbl),    pointer     :: swork(:)
    real(kind=dbl),    pointer     :: pmm(:), pmj2(:), pmj1(:), pmj(:)
    real(kind=dbl),    pointer     :: cosx(:), sinx(:), cosx2(:)
    real(kind=dbl),    pointer     :: sumN(:), sumS(:)
    complex(kind=dbl), allocatable :: rcc(:)
    type(c_ptr)                    :: c_swork, c_pmm, c_pmj2, c_pmj1, c_pmj, c_cosx, c_sinx, c_cosx2, c_sumN, c_sumS
    
    !Transform to suitable real input
    call this%lgp%alloc_cscal_sub( 1, rcc )
    call this%lgp%index_bwd_sub( 1, cc, rcc )
    
    !Allocating memory
    call calloc_sub( 4*step,                  c_swork ); call c_f_pointer( c_swork, swork, [ 4*step                  ] )
    call calloc_sub(   step,                  c_pmm   ); call c_f_pointer( c_pmm,   pmm,   [   step                  ] )
    call calloc_sub(   step,                  c_pmj   ); call c_f_pointer( c_pmj,   pmj,   [   step                  ] )
    call calloc_sub(   step,                  c_pmj1  ); call c_f_pointer( c_pmj1,  pmj1,  [   step                  ] )
    call calloc_sub(   step,                  c_pmj2  ); call c_f_pointer( c_pmj2,  pmj2,  [   step                  ] )
    call calloc_sub(   step,                  c_cosx  ); call c_f_pointer( c_cosx,  cosx,  [   step                  ] )
    call calloc_sub(   step,                  c_cosx2 ); call c_f_pointer( c_cosx2, cosx2, [   step                  ] )
    call calloc_sub(   step,                  c_sinx  ); call c_f_pointer( c_sinx,  sinx,  [   step                  ] )
    call calloc_sub(   step*this%fourtrans%n, c_sumN  ); call c_f_pointer( c_sumN,  sumN,  [   step*this%fourtrans%n ] )
    call calloc_sub(   step*this%fourtrans%n, c_sumS  ); call c_f_pointer( c_sumS,  sumS,  [   step*this%fourtrans%n ] )
    
    !Cycle over latitudes :: calculating step at once
    do itheta = 1, (this%lgp%nLege/step)*step, step
      cosx  = this%lgp%rw(itheta:itheta+step-1,1)
      sinx  = this%lgp%rw(itheta:itheta+step-1,2)
      cosx2 = this%lgp%rw(itheta:itheta+step-1,3)
      
      call zero_rarray_sub( step*this%fourtrans%n, sumN )
      call zero_rarray_sub( step*this%fourtrans%n, sumS )
      
      call this%lgp%bwd_legesum_sub( 1, rcc, sumN, sumS, cosx, sinx, cosx2, pmm, pmj2, pmj1, pmj, swork )
      
      call this%fourtrans%fft_c2r_sub( step, sumN )
      call this%fourtrans%fft_c2r_sub( step, sumS )
      
      do concurrent ( i2 = 1:this%fourtrans%n, i1 = 0:step-1 )
        grid(itheta+i1,i2,1) = sumN(i1+1+step*(i2-1))
        grid(itheta+i1,i2,2) = sumS(i1+1+step*(i2-1))
      end do
    end do
    
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
    
    deallocate( rcc )
    
  end procedure space_to_grid_sub
  
  module procedure grid_to_space_sub
    integer                        :: itheta, i1, i2
    real(kind=dbl),    pointer     :: swork(:)
    real(kind=dbl),    pointer     :: pmm(:), pmj2(:), pmj1(:), pmj(:)
    real(kind=dbl),    pointer     :: cosx(:), sinx(:), cosx2(:), wght(:)
    real(kind=dbl),    pointer     :: sumN(:), sumS(:)
    complex(kind=dbl), allocatable :: crr(:), rcr(:)
    type(c_ptr)                    :: c_swork, c_pmm, c_pmj2, c_pmj1, c_pmj, c_cosx, c_sinx, c_cosx2, c_wght, c_sumN, c_sumS
    
    !Allocate input array
    call this%lgp%alloc_cscal_sub( 1, rcr )
    call this%reindexing%allocate_scalars_sub( 1, crr )
    
    !Allocating memory
    call calloc_sub( 4*step,                  c_swork ); call c_f_pointer( c_swork, swork, [ 4*step                  ] )
    call calloc_sub(   step,                  c_pmm   ); call c_f_pointer( c_pmm,   pmm,   [   step                  ] )
    call calloc_sub(   step,                  c_pmj   ); call c_f_pointer( c_pmj,   pmj,   [   step                  ] )
    call calloc_sub(   step,                  c_pmj1  ); call c_f_pointer( c_pmj1,  pmj1,  [   step                  ] )
    call calloc_sub(   step,                  c_pmj2  ); call c_f_pointer( c_pmj2,  pmj2,  [   step                  ] )
    call calloc_sub(   step,                  c_cosx  ); call c_f_pointer( c_cosx,  cosx,  [   step                  ] )
    call calloc_sub(   step,                  c_cosx2 ); call c_f_pointer( c_cosx2, cosx2, [   step                  ] )
    call calloc_sub(   step,                  c_sinx  ); call c_f_pointer( c_sinx,  sinx,  [   step                  ] )
    call calloc_sub(   step,                  c_wght  ); call c_f_pointer( c_wght,  wght,  [   step                  ] )
    call calloc_sub(   step*this%fourtrans%n, c_sumN  ); call c_f_pointer( c_sumN,  sumN,  [   step*this%fourtrans%n ] )
    call calloc_sub(   step*this%fourtrans%n, c_sumS  ); call c_f_pointer( c_sumS,  sumS,  [   step*this%fourtrans%n ] )
    
    !Cycle over latitudes :: computing step at once
    do itheta = 1, (this%lgp%nLege/step)*step, step
      cosx  = this%lgp%rw(itheta:itheta+step-1,1)
      sinx  = this%lgp%rw(itheta:itheta+step-1,2)
      cosx2 = this%lgp%rw(itheta:itheta+step-1,3)
      wght  = this%lgp%rw(itheta:itheta+step-1,4)
      
      do concurrent ( i2 = 1:this%fourtrans%n, i1 = 0:step-1 )
        sumN(i1+1+step*(i2-1)) = grid(itheta+i1,i2,1)
        sumS(i1+1+step*(i2-1)) = grid(itheta+i1,i2,2)
      end do
      
      call this%fourtrans%fft_r2c_sub( step, sumN )
      call this%fourtrans%fft_r2c_sub( step, sumS )
      
      call this%lgp%fwd_legesum_sub( 1, sumN, sumS, rcr, cosx, sinx, cosx2, wght, pmm, pmj2, pmj1, pmj, swork )
    end do
    
    !Reindex input array
    call this%lgp%index_fwd_sub( 1, crr, rcr )
    call this%reindexing%scal2scal_mj_to_jm_sub( crr, 1, 1, cr, 1, 1 )
    
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
    
    deallocate( crr, rcr )
    
  end procedure grid_to_space_sub
  
end submodule sgs