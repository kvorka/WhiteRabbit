submodule (lateral_grid) transform
  implicit none; contains
  
  module procedure transform_sub
    integer                     :: itheta
    real(kind=dbl), allocatable :: sumN(:), sumS(:), rwork(:)
    
    !Allocating memory :: step is set to 16
    allocate( sumN(16*nb*this%fourtrans%n), sumS(16*nb*this%fourtrans%n), rwork(48) )
    
    !Cycle over latitudes :: calculating 16 at once
    do itheta = 1, (this%lgp%nLege/16)*16, 16
      rwork( 1:16) = this%lgp%rw(itheta:itheta+15,1)
      rwork(17:32) = this%lgp%rw(itheta:itheta+15,2)
      rwork(33:48) = this%lgp%rw(itheta:itheta+15,3)
      
      call zero_rarray_sub( 16*nb*this%fourtrans%n, sumN(1) )
      call zero_rarray_sub( 16*nb*this%fourtrans%n, sumS(1) )
      
      call this%lgp%backward_legesum_16_sub( nb, cc(1,1), sumN(1), sumS(1), rwork(1) )
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 16*nb, sumN(1) )
      call grid_sub( 16, this%fourtrans%n, sumN(1) )
      call this%fourtrans%fft_r2c_sub( 16*nf, sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 16*nb, sumS(1) )
      call grid_sub( 16, this%fourtrans%n, sumS(1) )
      call this%fourtrans%fft_r2c_sub( 16*nf, sumS(1) )
      
      call this%lgp%forward_legesum_16_sub( nf, sumN(1), sumS(1), cr(1,1), rwork(1) )
    end do
    
    !Cycle over latitudes :: calculating 8 at once
    do itheta = (this%lgp%nLege/16)*16+1, (this%lgp%nLege/8)*8, 8
      rwork( 1: 8) = this%lgp%rw(itheta:itheta+7,1)
      rwork( 9:16) = this%lgp%rw(itheta:itheta+7,2)
      rwork(17:24) = this%lgp%rw(itheta:itheta+7,3)
      
      call zero_rarray_sub( 8*nb*this%fourtrans%n, sumN(1) )
      call zero_rarray_sub( 8*nb*this%fourtrans%n, sumS(1) )
      
      call this%lgp%backward_legesum_8_sub( nb, cc(1,1), sumN(1), sumS(1), rwork(1) )
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 8*nb, sumN(1) )
      call grid_sub( 8, this%fourtrans%n, sumN(1) )
      call this%fourtrans%fft_r2c_sub( 8*nf, sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 8*nb, sumS(1) )
      call grid_sub( 8, this%fourtrans%n, sumS(1) )
      call this%fourtrans%fft_r2c_sub( 8*nf, sumS(1) )
      
      call this%lgp%forward_legesum_8_sub( nf, sumN(1), sumS(1), cr(1,1), rwork(1) )
    end do
    
    !Cycle over latitudes :: calculating 4 possibly left at once
    do itheta = (this%lgp%nLege/8)*8+1, this%lgp%nLege, 4
      rwork(1: 4) = this%lgp%rw(itheta:itheta+3,1)
      rwork(5: 8) = this%lgp%rw(itheta:itheta+3,2)
      rwork(9:12) = this%lgp%rw(itheta:itheta+3,3)
      
      call zero_rarray_sub( 4*nb*this%fourtrans%n, sumN(1) )
      call zero_rarray_sub( 4*nb*this%fourtrans%n, sumS(1) )
      
      !Sums of associate Legendre polynomials
      call this%lgp%backward_legesum_4_sub( nb, cc(1,1), sumN(1), sumS(1), rwork(1) )
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 4*nb, sumN(1) )
      call grid_sub( 4, this%fourtrans%n, sumN(1) )
      call this%fourtrans%fft_r2c_sub( 4*nf, sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 4*nb, sumS(1) )
      call grid_sub( 4, this%fourtrans%n, sumS(1) )
      call this%fourtrans%fft_r2c_sub( 4*nf, sumS(1) )
      
      !Weighted sums of associate Legendre polynomials
      call this%lgp%forward_legesum_4_sub( nf, sumN(1), sumS(1), cr(1,1), rwork(1) )
    end do
    
    !Cleaning
    deallocate( sumN, sumS, rwork )
    
  end procedure transform_sub
  
end submodule transform