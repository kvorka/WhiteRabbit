submodule (lateral_grid) sgs
  implicit none; contains
  
  module procedure allocate_grid_sub
    
    if ( allocated(grid) ) deallocate( grid )
    
    allocate( grid(this%lgp%nLege,this%fourtrans%n,2) )
      call zero_rarray_sub( 2 * this%lgp%nLege * this%fourtrans%n, grid )
    
  end procedure allocate_grid_sub
  
  module procedure space_to_grid_sub
    integer                     :: itheta, i1, i2
    real(kind=dbl), allocatable :: sumN(:), sumS(:), rwork(:)
    
    !Allocating memory
    allocate( sumN(16*this%fourtrans%n), sumS(16*this%fourtrans%n), rwork(32) )
    
    !Cycle over latitudes :: calculating 16 at once
    do itheta = 1, (this%lgp%nLege/16)*16, 16
      rwork( 1:16) = this%lgp%rw(itheta:itheta+15,1)
      rwork(17:32) = this%lgp%rw(itheta:itheta+15,2)
      
      call zero_rarray_sub( 16*this%fourtrans%n, sumN(1) )
      call zero_rarray_sub( 16*this%fourtrans%n, sumS(1) )
      
      call this%lgp%backward_legesum_16_sub( 1, cc(1), sumN(1), sumS(1), rwork(1) )
      
      call this%fourtrans%fft_c2r_sub( 16, sumN(1) )
      call this%fourtrans%fft_c2r_sub( 16, sumS(1) )
      
      do concurrent ( i2 = 1:this%fourtrans%n, i1 = 0:15 )
        grid(itheta+i1,i2,1) = sumN(i1+1+16*(i2-1))
        grid(itheta+i1,i2,2) = sumS(i1+1+16*(i2-1))
      end do
    end do
    
    !Cycle over latitudes :: calculating 8 at once
    do itheta = (this%lgp%nLege/16)*16+1, (this%lgp%nLege/8)*8, 8
      rwork(1: 8) = this%lgp%rw(itheta:itheta+7,1)
      rwork(9:16) = this%lgp%rw(itheta:itheta+7,2)
      
      call zero_rarray_sub( 8*this%fourtrans%n, sumN(1) )
      call zero_rarray_sub( 8*this%fourtrans%n, sumS(1) )
      
      call this%lgp%backward_legesum_8_sub( 1, cc(1), sumN(1), sumS(1), rwork(1) )
      
      call this%fourtrans%fft_c2r_sub( 8, sumN(1) )
      call this%fourtrans%fft_c2r_sub( 8, sumS(1) )
      
      do concurrent ( i2 = 1:this%fourtrans%n, i1 = 0:7 )
        grid(itheta+i1,i2,1) = sumN(i1+1+8*(i2-1))
        grid(itheta+i1,i2,2) = sumS(i1+1+8*(i2-1))
      end do
    end do
    
    !Cycle over latitudes :: calculating 4 at once
    do itheta = (this%lgp%nLege/8)*8+1, this%lgp%nLege, 4
      rwork(1:4) = this%lgp%rw(itheta:itheta+3,1)
      rwork(5:8) = this%lgp%rw(itheta:itheta+3,2)
      
      call zero_rarray_sub( 4*this%fourtrans%n, sumN(1) )
      call zero_rarray_sub( 4*this%fourtrans%n, sumS(1) )
      
      call this%lgp%backward_legesum_4_sub( 1, cc(1), sumN(1), sumS(1), rwork(1) )
      
      call this%fourtrans%fft_c2r_sub( 4, sumN(1) )
      call this%fourtrans%fft_c2r_sub( 4, sumS(1) )
      
      do concurrent ( i2 = 1:this%fourtrans%n, i1 = 0:3 )
        grid(itheta+i1,i2,1) = sumN(i1+1+4*(i2-1))
        grid(itheta+i1,i2,2) = sumS(i1+1+4*(i2-1))
      end do
    end do
    
    !Cleaning
    deallocate( sumN, sumS, rwork )
    
  end procedure space_to_grid_sub
  
  module procedure grid_to_space_sub
    integer                     :: itheta, i1, i2
    real(kind=dbl), allocatable :: sumN(:), sumS(:), rwork(:)
    
    !Allocating memory
    allocate( sumN(16*this%fourtrans%n), sumS(16*this%fourtrans%n), rwork(48) )
    
    !Cycle over latitudes :: computing 16 at once
    do itheta = 1, (this%lgp%nLege/16)*16, 16
      rwork( 1:16) = this%lgp%rw(itheta:itheta+15,1)
      rwork(17:32) = this%lgp%rw(itheta:itheta+15,2)
      rwork(33:48) = this%lgp%rw(itheta:itheta+15,3)
      
      do concurrent ( i2 = 1:this%fourtrans%n, i1 = 0:15 )
        sumN(i1+1+16*(i2-1)) = grid(itheta+i1,i2,1)
        sumS(i1+1+16*(i2-1)) = grid(itheta+i1,i2,2)
      end do
      
      call this%fourtrans%fft_r2c_sub( 16, sumN(1) )
      call this%fourtrans%fft_r2c_sub( 16, sumS(1) )
      
      call this%lgp%forward_legesum_16_sub( 1, sumN(1), sumS(1), cr(1), rwork(1) )
    end do
    
    !Cycle over latitudes :: computing 8 at once
    do itheta = (this%lgp%nLege/16)*16+1, (this%lgp%nLege/8)*8, 8
      rwork( 1: 8) = this%lgp%rw(itheta:itheta+7,1)
      rwork( 9:16) = this%lgp%rw(itheta:itheta+7,2)
      rwork(17:24) = this%lgp%rw(itheta:itheta+7,3)
      
      do concurrent ( i2 = 1:this%fourtrans%n, i1 = 0:7 )
        sumN(i1+1+8*(i2-1)) = grid(itheta+i1,i2,1)
        sumS(i1+1+8*(i2-1)) = grid(itheta+i1,i2,2)
      end do
      
      call this%fourtrans%fft_r2c_sub( 8, sumN(1) )
      call this%fourtrans%fft_r2c_sub( 8, sumS(1) )
      
      call this%lgp%forward_legesum_8_sub( 1, sumN(1), sumS(1), cr(1), rwork(1) )
    end do
    
    !Cycle over latitudes :: computing 4 at once
    do itheta = (this%lgp%nLege/8)*8+1, this%lgp%nLege, 4
      rwork(1: 4) = this%lgp%rw(itheta:itheta+3,1)
      rwork(5: 8) = this%lgp%rw(itheta:itheta+3,2)
      rwork(9:12) = this%lgp%rw(itheta:itheta+3,3)
      
      do concurrent ( i2 = 1:this%fourtrans%n, i1 = 0:3 )
        sumN(i1+1+4*(i2-1)) = grid(itheta+i1,i2,1)
        sumS(i1+1+4*(i2-1)) = grid(itheta+i1,i2,2)
      end do
      
      call this%fourtrans%fft_r2c_sub( 4, sumN(1) )
      call this%fourtrans%fft_r2c_sub( 4, sumS(1) )
      
      call this%lgp%forward_legesum_4_sub( 1, sumN(1), sumS(1), cr(1), rwork(1) )
    end do
    
    !Cleaning
    deallocate( sumN, sumS, rwork )
    
  end procedure grid_to_space_sub
  
end submodule sgs