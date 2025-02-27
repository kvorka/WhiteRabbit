submodule (lateral_grid) sgs
  implicit none; contains
  
  module procedure allocate_grid_sub
    
    if ( allocated(grid) ) deallocate( grid )
    
    allocate( grid(this%lgp%nLege,this%fourtrans%n,2) )
      call zero_rarray_sub( 2 * this%lgp%nLege * this%fourtrans%n, grid )
    
  end procedure allocate_grid_sub
  
  module procedure space_to_grid_sub
    integer                             :: itheta
    type(c_ptr)                         :: c_work, c_rcc
    real(kind=dbl), contiguous, pointer :: work(:)
    real(kind=dbl), contiguous, pointer :: rcc(:)
    real(kind=dbl), contiguous, pointer :: pwork(:)
    real(kind=dbl), contiguous, pointer :: cosx(:), sinx(:), cosx2(:)
    real(kind=dbl), contiguous, pointer :: sumN(:), sumS(:), swork(:)
    
    !Transform to suitable real input
    call this%lgp%alloc_rscal_sub( 1, c_rcc, rcc )
    call this%lgp%index_bwd_sub( 1, cc, rcc )
    
    !Allocating memory
    call alloc_aligned1d_sub( (7+2*this%fourtrans%n)*step, c_work, work )
      
    pwork => work(                             1 : ( 3                    ) * step )
    swork => work(   3*                   step+1 : ( 7                    ) * step )
    sumN  => work(   7*                   step+1 : ( 7+  this%fourtrans%n ) * step )
    sumS  => work( ( 7+this%fourtrans%n )*step+1 : ( 7+2*this%fourtrans%n ) * step )
    
    !Cycle over latitudes :: calculating step at once
    do itheta = 1, (this%lgp%nLege/step)*step, step
      cosx  => this%lgp%rw(itheta:itheta+step-1,1)
      sinx  => this%lgp%rw(itheta:itheta+step-1,2)
      cosx2 => this%lgp%rw(itheta:itheta+step-1,3)
      
      call zero_rarray_sub( step*this%fourtrans%n, sumN )
      call zero_rarray_sub( step*this%fourtrans%n, sumS )
      
      call this%lgp%bwd_legesum_sub( 1, rcc, sumN, sumS, cosx, sinx, cosx2, pwork, swork )
      
      call this%fourtrans%fft_c2r_sub( step, sumN )
      call this%fourtrans%fft_c2r_sub( step, sumS )
      
      call grid_op_save_sub( this%fourtrans%n, sumN, grid(itheta,1,1), this%lgp%nLege )
      call grid_op_save_sub( this%fourtrans%n, sumS, grid(itheta,1,2), this%lgp%nLege )
    end do
    
    !Cleaning
    call free_aligned1d_sub( c_work, work )
    call free_aligned1d_sub( c_rcc, rcc )
    
  end procedure space_to_grid_sub
  
  module procedure grid_to_space_sub
    integer                             :: itheta
    type(c_ptr)                         :: c_work, c_rcr
    real(kind=dbl), contiguous, pointer :: work(:)
    real(kind=dbl), contiguous, pointer :: rcr(:)
    real(kind=dbl), contiguous, pointer :: pwork(:)
    real(kind=dbl), contiguous, pointer :: cosx(:), sinx(:), cosx2(:), wght(:)
    real(kind=dbl), contiguous, pointer :: sumN(:), sumS(:), swork(:)
    complex(kind=dbl), allocatable      :: crr(:)
    
    !Allocate input array
    call this%lgp%alloc_rscal_sub( 1, c_rcr, rcr )
    call this%reindexing%allocate_scalars_sub( 1, crr )
    
    !Allocating memory
    call alloc_aligned1d_sub( (7+2*this%fourtrans%n)*step, c_work, work )
    
    pwork => work(                             1 : ( 3                    ) * step )
    swork => work(   3*                   step+1 : ( 7                    ) * step )
    sumN  => work(   7*                   step+1 : ( 7+  this%fourtrans%n ) * step )
    sumS  => work( ( 7+this%fourtrans%n )*step+1 : ( 7+2*this%fourtrans%n ) * step )
    
    !Cycle over latitudes :: computing step at once
    do itheta = 1, (this%lgp%nLege/step)*step, step
      cosx  => this%lgp%rw(itheta:itheta+step-1,1)
      sinx  => this%lgp%rw(itheta:itheta+step-1,2)
      cosx2 => this%lgp%rw(itheta:itheta+step-1,3)
      wght  => this%lgp%rw(itheta:itheta+step-1,4)
      
      call grid_op_load_sub( this%fourtrans%n, sumN, grid(itheta,1,1), this%lgp%nLege )
      call grid_op_load_sub( this%fourtrans%n, sumS, grid(itheta,1,2), this%lgp%nLege )
      
      call this%fourtrans%fft_r2c_sub( step, sumN )
      call this%fourtrans%fft_r2c_sub( step, sumS )
      
      call this%lgp%fwd_legesum_sub( 1, sumN, sumS, rcr, cosx, sinx, cosx2, wght, pwork, swork )
    end do
    
    !Reindex input array
    call this%lgp%index_fwd_sub( 1, crr, rcr )
    call this%reindexing%scal2scal_mj_to_jm_sub( crr, 1, 1, cr, 1, 1 )
    
    !Cleaning
    call free_aligned1d_sub( c_work, work )
    call free_aligned1d_sub( c_rcr, rcr   )
    
    deallocate( crr )
    
  end procedure grid_to_space_sub
  
end submodule sgs