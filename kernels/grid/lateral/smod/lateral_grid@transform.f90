submodule (lateral_grid) transform
  implicit none; contains
  
  module pure subroutine transform_sub(this, nf, nb, cc, cr, grid_sub)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: nf, nb
    complex(kind=dbl),    intent(in)    :: cc(nb,*)
    complex(kind=dbl),    intent(inout) :: cr(nf,*)
    integer                             :: i
    real(kind=dbl),       allocatable   :: grid(:)
    complex(kind=dbl),    allocatable   :: sumN(:), sumS(:)
    
    interface
      module pure subroutine grid_sub(step, nfour, gxyz)
        integer,                intent(in)    :: step, nfour
        real(kind=dbl), target, intent(inout) :: gxyz(*)
      end subroutine grid_sub
    end interface
    
    !Allocating memory :: step is set to 16
    allocate( sumN(16*nb*(this%lgp%jmax+1)), sumS(16*nb*(this%lgp%jmax+1)), grid(16*nb*this%nFourier) )
    
    !Cycle over latitudes :: calculating 16 at once
    do i = 1, (this%nLegendre/16)*16, 16
      call this%lgp%backward_legesum_16_sub( nb, this%cosx(i), cc(1,1), sumN(1), sumS(1) )
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 16*nb, sumN(1), grid(1) )
      call grid_sub( 16, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 16*nf, grid(1), sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 16*nb, sumS(1), grid(1) )
      call grid_sub( 16, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 16*nf, grid(1), sumS(1) )
      
      call this%lgp%forward_legesum_16_sub( nf, this%cosx(i), this%weight(i), sumN(1), sumS(1), cr(1,1) )
    end do

    !Cycle over latitudes :: calculating 8 at once
    do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
      call this%lgp%backward_legesum_8_sub( nb, this%cosx(i), cc(1,1), sumN(1), sumS(1) )
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 8*nb, sumN(1), grid(1) )
      call grid_sub( 8, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 8*nf, grid(1), sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 8*nb, sumS(1), grid(1) )
      call grid_sub( 8, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 8*nf, grid(1), sumS(1) )
      
      call this%lgp%forward_legesum_8_sub( nf, this%cosx(i), this%weight(i), sumN(1), sumS(1), cr(1,1) )
    end do
    
    !Cycle over latitudes :: calculating 4 possibly left at once
    do i = (this%nLegendre/8)*8+1, this%nLegendre, 4
      !Sums of associate Legendre polynomials
      call this%lgp%backward_legesum_4_sub( nb, this%cosx(i), cc(1,1), sumN(1), sumS(1) )
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 4*nb, sumN(1), grid(1) )
      call grid_sub( 4, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 4*nf, grid(1), sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 4*nb, sumS(1), grid(1) )
      call grid_sub( 4, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 4*nf, grid(1), sumS(1) )
      
      !Weighted sums of associate Legendre polynomials
      call this%lgp%forward_legesum_4_sub( nf, this%cosx(i), this%weight(i), sumN(1), sumS(1), cr(1,1) )
    end do
    
    !Cleaning
    deallocate( sumN, sumS, grid )
    
  end subroutine transform_sub
  
end submodule transform