submodule (lateral_grid) sgs
  implicit none; contains
  
  module pure subroutine space_to_grid_sub(this, cc, grid)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cc(*)
    real(kind=dbl),       intent(out) :: grid(:,:,:)
    integer                           :: itheta
    complex(kind=dbl),    allocatable :: sumN(:), sumS(:)
    
    !Allocating memory
    allocate( sumN(16*(this%lgp%jmax+1)), sumS(16*(this%lgp%jmax+1)) )
    
    !Cycle over latitudes :: calculating 16 at once
    do itheta = 1, (this%lgp%nLege/16)*16, 16
      call this%lgp%backward_legesum_16_sub( itheta, 1, cc(1), sumN(1), sumS(1) )
      
      call this%fourtrans%exec_c2r_sub( 16, sumN(1), grid(itheta:itheta+15,:,1) )
      call this%fourtrans%exec_c2r_sub( 16, sumS(1), grid(itheta:itheta+15,:,2) )
    end do
    
    !Cycle over latitudes :: calculating 8 at once
    do itheta = (this%lgp%nLege/16)*16+1, (this%lgp%nLege/8)*8, 8
      call this%lgp%backward_legesum_8_sub( itheta, 1, cc(1), sumN(1), sumS(1) )
      
      call this%fourtrans%exec_c2r_sub( 8, sumN(1), grid(itheta:itheta+7,:,1) )
      call this%fourtrans%exec_c2r_sub( 8, sumS(1), grid(itheta:itheta+7,:,2) )
    end do
    
    !Cycle over latitudes :: calculating 4 at once
    do itheta = (this%lgp%nLege/8)*8+1, this%lgp%nLege, 4
      call this%lgp%backward_legesum_4_sub( itheta, 1, cc(1), sumN(1), sumS(1) )
      
      call this%fourtrans%exec_c2r_sub( 4, sumN(1), grid(itheta:itheta+3,:,1) )
      call this%fourtrans%exec_c2r_sub( 4, sumS(1), grid(itheta:itheta+3,:,2) )
    end do
    
    !Cleaning
    deallocate( sumN, sumS )
    
  end subroutine space_to_grid_sub
  
  module pure subroutine grid_to_space_sub(this, grid, cr)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(inout) :: grid(:,:,:)
    complex(kind=dbl),    intent(out)   :: cr(*)
    integer                             :: itheta
    complex(kind=dbl),    allocatable   :: sumN(:), sumS(:)
    
    !Allocating memory
    allocate( sumN(16*(this%lgp%jmax+1)), sumS(16*(this%lgp%jmax+1)) )
    
    !Cycle over latitudes :: computing 16 at once
    do itheta = 1, (this%lgp%nLege/16)*16, 16
      call this%fourtrans%exec_r2c_sub( 16, grid(itheta:itheta+15,:,1), sumN(1) )
      call this%fourtrans%exec_r2c_sub( 16, grid(itheta:itheta+15,:,2), sumS(1) )
      
      call this%lgp%forward_legesum_16_sub( itheta, 1, sumN(1), sumS(1), cr(1) )
    end do
    
    !Cycle over latitudes :: computing 8 at once
    do itheta = (this%lgp%nLege/16)*16+1, (this%lgp%nLege/8)*8, 8
      call this%fourtrans%exec_r2c_sub( 8, grid(itheta:itheta+7,:,1), sumN(1) )
      call this%fourtrans%exec_r2c_sub( 8, grid(itheta:itheta+7,:,2), sumS(1) )
      
      call this%lgp%forward_legesum_8_sub( itheta, 1, sumN(1), sumS(1), cr(1) )
    end do
    
    !Cycle over latitudes :: computing 4 at once
    do itheta = (this%lgp%nLege/8)*8+1, this%lgp%nLege, 4
      call this%fourtrans%exec_r2c_sub( 4, grid(itheta:itheta+3,:,1), sumN(1) )
      call this%fourtrans%exec_r2c_sub( 4, grid(itheta:itheta+3,:,2), sumS(1) )
      
      call this%lgp%forward_legesum_4_sub( itheta, 1, sumN(1), sumS(1), cr(1) )
    end do
    
    !Cleaning
    deallocate( sumN, sumS )
    
  end subroutine grid_to_space_sub
  
end submodule sgs