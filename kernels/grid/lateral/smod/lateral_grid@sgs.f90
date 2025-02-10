submodule (lateral_grid) sgs
  implicit none; contains
  
  module procedure allocate_grid_sub
    
    if ( allocated(grid) ) deallocate( grid )
    
    allocate( grid(this%lgp%nLege,this%fourtrans%n,2) )
      call zero_rarray_sub( 2 * this%lgp%nLege * this%fourtrans%n, grid )
    
  end procedure allocate_grid_sub
  
  module procedure space_to_grid_sub
    integer                        :: itheta, i1, i2
    real(kind=dbl),    pointer     :: swork(:), pmm(:), pmj2(:), pmj1(:), pmj(:)
    real(kind=dbl),    pointer     :: cosx(:), sinx(:), cosx2(:)
    real(kind=dbl),    pointer     :: sumN(:), sumS(:)
    complex(kind=dbl), allocatable :: rcc(:)
    type(c_ptr)                    :: c_swork, c_pmm, c_pmj2, c_pmj1, c_pmj
    type(c_ptr)                    :: c_cosx, c_sinx, c_cosx2
    type(c_ptr)                    :: c_sumN, c_sumS
    
    !Transform to suitable real input
    call this%lgp%alloc_cscal_sub( 1, rcc )
    call this%lgp%index_bwd_sub( 1, cc, rcc )
    
    !Allocating memory
    call alloc_aligned_sub( 4*step,                  c_swork, swork )
    call alloc_aligned_sub(   step,                  c_pmm,   pmm   )
    call alloc_aligned_sub(   step,                  c_pmj,   pmj   )
    call alloc_aligned_sub(   step,                  c_pmj1,  pmj1  )
    call alloc_aligned_sub(   step,                  c_pmj2,  pmj2  )
    call alloc_aligned_sub(   step,                  c_cosx,  cosx  )
    call alloc_aligned_sub(   step,                  c_cosx2, cosx2 )
    call alloc_aligned_sub(   step,                  c_sinx,  sinx  )
    call alloc_aligned_sub(   step*this%fourtrans%n, c_sumN,  sumN  )
    call alloc_aligned_sub(   step*this%fourtrans%n, c_sumS,  sumS  )
    
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
    call free_aligned_sub( c_swork, swork )
    call free_aligned_sub( c_pmm,   pmm   )
    call free_aligned_sub( c_pmj,   pmj   )
    call free_aligned_sub( c_pmj1,  pmj1  )
    call free_aligned_sub( c_pmj2,  pmj2  )
    call free_aligned_sub( c_cosx,  cosx  )
    call free_aligned_sub( c_sinx,  sinx  )
    call free_aligned_sub( c_cosx2, cosx2 )
    call free_aligned_sub( c_sumN,  sumN  )
    call free_aligned_sub( c_sumS,  sumS  )
    
    deallocate( rcc )
    
  end procedure space_to_grid_sub
  
  module procedure grid_to_space_sub
    integer                        :: itheta, i1, i2
    real(kind=dbl),    pointer     :: swork(:), pmm(:), pmj2(:), pmj1(:), pmj(:)
    real(kind=dbl),    pointer     :: cosx(:), sinx(:), cosx2(:), wght(:)
    real(kind=dbl),    pointer     :: sumN(:), sumS(:)
    complex(kind=dbl), allocatable :: crr(:), rcr(:)
    type(c_ptr)                    :: c_swork, c_pmm, c_pmj2, c_pmj1, c_pmj
    type(c_ptr)                    :: c_cosx, c_sinx, c_cosx2, c_wght
    type(c_ptr)                    :: c_sumN, c_sumS
    
    !Allocate input array
    call this%lgp%alloc_cscal_sub( 1, rcr )
    call this%reindexing%allocate_scalars_sub( 1, crr )
    
    !Allocating memory
    call alloc_aligned_sub( 4*step,                  c_swork, swork )
    call alloc_aligned_sub(   step,                  c_pmm,   pmm   )
    call alloc_aligned_sub(   step,                  c_pmj,   pmj   )
    call alloc_aligned_sub(   step,                  c_pmj1,  pmj1  )
    call alloc_aligned_sub(   step,                  c_pmj2,  pmj2  )
    call alloc_aligned_sub(   step,                  c_cosx,  cosx  )
    call alloc_aligned_sub(   step,                  c_cosx2, cosx2 )
    call alloc_aligned_sub(   step,                  c_sinx,  sinx  )
    call alloc_aligned_sub(   step,                  c_wght,  wght  )
    call alloc_aligned_sub(   step*this%fourtrans%n, c_sumN,  sumN  )
    call alloc_aligned_sub(   step*this%fourtrans%n, c_sumS,  sumS  )
    
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
    
    deallocate( crr, rcr )
    
  end procedure grid_to_space_sub
  
end submodule sgs