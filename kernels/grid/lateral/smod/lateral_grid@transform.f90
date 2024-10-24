submodule (lateral_grid) transform
  implicit none; contains
  
  module pure subroutine transform_sub(this, nf, nb, cc, cr, grid_sub)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: nf, nb
    complex(kind=dbl),    intent(in)    :: cc(nb,*)
    complex(kind=dbl),    intent(inout) :: cr(nf,*)
    integer                             :: itheta
    real(kind=dbl),       allocatable   :: sumN(:), sumS(:)
    
    interface
      pure subroutine grid_sub(step, nfour, gxyz); import :: dbl
        integer,                intent(in)    :: step, nfour
        real(kind=dbl), target, intent(inout) :: gxyz(step,*)
      end subroutine grid_sub
    end interface
    
    !Allocating memory :: step is set to 16
    allocate( sumN(16*nb*this%fourtrans%n), sumS(16*nb*this%fourtrans%n) )
    
    !Cycle over latitudes :: calculating 16 at once
    do itheta = 1, (this%lgp%nLege/16)*16, 16
      call zero_rarray_sub(16*nb*this%fourtrans%n, sumN(1))
      call zero_rarray_sub(16*nb*this%fourtrans%n, sumS(1))
      
      call this%lgp%backward_legesum_16_sub( itheta, nb, cc(1,1), sumN(1), sumS(1) )
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 16*nb, sumN(1) )
      call grid_sub( 16, this%fourtrans%n, sumN(1) )
      call this%fourtrans%fft_r2c_sub( 16*nf, sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 16*nb, sumS(1) )
      call grid_sub( 16, this%fourtrans%n, sumS(1) )
      call this%fourtrans%fft_r2c_sub( 16*nf, sumS(1) )
      
      call this%lgp%forward_legesum_16_sub( itheta, nf, sumN(1), sumS(1), cr(1,1) )
    end do
    
    !Cycle over latitudes :: calculating 8 at once
    do itheta = (this%lgp%nLege/16)*16+1, (this%lgp%nLege/8)*8, 8
      call zero_rarray_sub(8*nb*this%fourtrans%n, sumN(1))
      call zero_rarray_sub(8*nb*this%fourtrans%n, sumS(1))
      
      call this%lgp%backward_legesum_8_sub( itheta, nb, cc(1,1), sumN(1), sumS(1) )
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 8*nb, sumN(1) )
      call grid_sub( 8, this%fourtrans%n, sumN(1) )
      call this%fourtrans%fft_r2c_sub( 8*nf, sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 8*nb, sumS(1) )
      call grid_sub( 8, this%fourtrans%n, sumS(1) )
      call this%fourtrans%fft_r2c_sub( 8*nf, sumS(1) )
      
      call this%lgp%forward_legesum_8_sub( itheta, nf, sumN(1), sumS(1), cr(1,1) )
    end do
    
    !Cycle over latitudes :: calculating 4 possibly left at once
    do itheta = (this%lgp%nLege/8)*8+1, this%lgp%nLege, 4
      call zero_rarray_sub(4*nb*this%fourtrans%n, sumN(1))
      call zero_rarray_sub(4*nb*this%fourtrans%n, sumS(1))
      
      !Sums of associate Legendre polynomials
      call this%lgp%backward_legesum_4_sub( itheta, nb, cc(1,1), sumN(1), sumS(1) )
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 4*nb, sumN(1))
      call grid_sub( 4, this%fourtrans%n, sumN(1) )
      call this%fourtrans%fft_r2c_sub( 4*nf, sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%fft_c2r_sub( 4*nb, sumS(1) )
      call grid_sub( 4, this%fourtrans%n, sumS(1) )
      call this%fourtrans%fft_r2c_sub( 4*nf, sumS(1) )
      
      !Weighted sums of associate Legendre polynomials
      call this%lgp%forward_legesum_4_sub( itheta, nf, sumN(1), sumS(1), cr(1,1) )
    end do
    
    !Cleaning
    deallocate( sumN, sumS )
    
  end subroutine transform_sub
  
end submodule transform