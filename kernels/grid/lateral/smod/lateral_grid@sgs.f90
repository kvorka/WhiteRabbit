submodule (lateral_grid) sgs
  implicit none; contains
  
  module procedure allocate_grid_sub
    
    if ( allocated(grid) ) deallocate( grid )
    
    allocate( grid(this%lgp%nLege,this%fourtrans%n,2) )
      call zero_rarray_sub( 2 * this%lgp%nLege * this%fourtrans%n, grid )
    
  end procedure allocate_grid_sub
  
  module procedure space_to_grid_sub
    integer                             :: itheta, i1, i2
    type(c_ptr)                         :: c_work
    real(kind=dbl), pointer, contiguous :: work(:)
    real(kind=dbl), pointer, contiguous :: pmm(:), pmj2(:), pmj1(:), pmj(:)
    real(kind=dbl), pointer, contiguous :: cosx(:), sinx(:), cosx2(:)
    real(kind=dbl), pointer, contiguous :: sumN(:), sumS(:), swork(:)
    real(kind=dbl), allocatable         :: rcc(:)
    
    !Transform to suitable real input
    call this%lgp%alloc_rscal_sub( 1, rcc )
    call this%lgp%index_bwd_sub( 1, cc, rcc )
    
    !Allocating memory
    call alloc_aligned1d_sub( 2*(4+this%fourtrans%n)*step, c_work, work )
      
      pmm   => work(                             1 :                          step )
      pmj   => work(                        step+1 :   2*                     step )
      pmj1  => work(   2*                   step+1 :   3*                     step )
      pmj2  => work(   3*                   step+1 :   4*                     step )
      swork => work(   4*                   step+1 :   8*                     step )
      sumN  => work(   8*                   step+1 : ( 8+  this%fourtrans%n )*step )
      sumS  => work( ( 8+this%fourtrans%n )*step+1 : ( 8+2*this%fourtrans%n )*step )
    
    !Cycle over latitudes :: calculating step at once
    do itheta = 1, (this%lgp%nLege/step)*step, step
      cosx  => this%lgp%rw(itheta:itheta+step-1,1)
      sinx  => this%lgp%rw(itheta:itheta+step-1,2)
      cosx2 => this%lgp%rw(itheta:itheta+step-1,3)
      
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
    call free_aligned1d_sub( c_work, work )
    
    deallocate( rcc )
    
  end procedure space_to_grid_sub
  
  module procedure grid_to_space_sub
    integer                                :: itheta, i1, i2
    type(c_ptr)                            :: c_work
    real(kind=dbl),    pointer, contiguous :: work(:)
    real(kind=dbl),    pointer, contiguous :: pmm(:), pmj2(:), pmj1(:), pmj(:)
    real(kind=dbl),    pointer, contiguous :: cosx(:), sinx(:), cosx2(:), wght(:)
    real(kind=dbl),    pointer, contiguous :: sumN(:), sumS(:), swork(:)
    real(kind=dbl),    allocatable         :: rcr(:)
    complex(kind=dbl), allocatable         :: crr(:)
    
    !Allocate input array
    call this%lgp%alloc_rscal_sub( 1, rcr )
    call this%reindexing%allocate_scalars_sub( 1, crr )
    
    !Allocating memory
    call alloc_aligned1d_sub( 2*(4+this%fourtrans%n)*step, c_work, work )
      
      pmm   => work(                             1 :                          step )
      pmj   => work(                        step+1 :   2*                     step )
      pmj1  => work(   2*                   step+1 :   3*                     step )
      pmj2  => work(   3*                   step+1 :   4*                     step )
      swork => work(   4*                   step+1 :   8*                     step )
      sumN  => work(   8*                   step+1 : ( 8+  this%fourtrans%n )*step )
      sumS  => work( ( 8+this%fourtrans%n )*step+1 : ( 8+2*this%fourtrans%n )*step )
    
    !Cycle over latitudes :: computing step at once
    do itheta = 1, (this%lgp%nLege/step)*step, step
      cosx  => this%lgp%rw(itheta:itheta+step-1,1)
      sinx  => this%lgp%rw(itheta:itheta+step-1,2)
      cosx2 => this%lgp%rw(itheta:itheta+step-1,3)
      wght  => this%lgp%rw(itheta:itheta+step-1,4)
      
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
    call free_aligned1d_sub( c_work, work )
    
    deallocate( crr, rcr )
    
  end procedure grid_to_space_sub
  
end submodule sgs