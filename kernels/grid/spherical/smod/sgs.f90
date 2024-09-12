submodule (SphericalHarmonics) sgs
  implicit none; contains
  
  module pure subroutine space_to_grid_sub(this, cc, grid)
    class(T_lateralGrid), intent(in)    :: this
    complex(kind=dbl),    intent(in)    :: cc(*)
    real(kind=dbl),       intent(out)   :: grid(:,:,:)
    integer                             :: i1, i2, i, j, m, mj
    real(kind=dbl),       allocatable   :: pmm(:), pmj0(:), pmj1(:), pmj2(:), csx(:), snx(:)
    complex(kind=dbl),    allocatable   :: ssm(:), asm(:), sumN(:), sumS(:)
    
    !Allocating memory :: step is set to 4
    allocate( pmm(4), pmj0(4), pmj1(4), pmj2(4), csx(4), snx(4), ssm(4), asm(4), &
            & sumN(4*(this%jmax+1)), sumS(4*(this%jmax+1))                       )
    
    !Cycle over latitudes :: step is set to 4
    do i = 1, this%nLegendre, 4
      do concurrent ( i2 = 1:4 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 4*(this%jmax+1), sumN(1) )
      call zero_carray_sub( 4*(this%jmax+1), sumS(1) )
      
      do m = 0, this%jmax
        call zero_carray_sub( 4, ssm(1) )
        call zero_carray_sub( 4, asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_4_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( 1, pmj0(1), cc(mj), ssm(1) )
          
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_4_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( 1, pmj0(1), cc(mj-1), asm(1) )
          
          call pmj_recursion_4_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( 1, pmj0(1), cc(mj), ssm(1) )
        end do
        
        if ( mod((this%jmax-m),2) /= 0 ) then
          call pmj_recursion_4_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( 1, pmj0(1), cc(mj+1), asm(1) )
        end if
        
        call pmj_backward_rcb_4_sub( 1, ssm(1), asm(1), sumN(4*m+1), sumS(4*m+1) )
      end do
      
      !Northern hemisphere :: fft towards grid
      call this%fourtrans%exec_c2r_sub( 4, sumN(1), grid(i:i+3,:,1) )
      
      !Southern hemisphere :: fft towards grid
      call this%fourtrans%exec_c2r_sub( 4, sumS(1), grid(i:i+3,:,2) )
    end do
    
    !Cleaning
    deallocate( pmm, pmj0, pmj1, pmj2, csx, snx, ssm, asm, sumN, sumS )
    
  end subroutine space_to_grid_sub
  
  module pure subroutine grid_to_space_sub(this, grid, cr)
    class(T_lateralGrid), intent(in)    :: this
    real(kind=dbl),       intent(inout) :: grid(:,:,:)
    complex(kind=dbl),    intent(out)   :: cr(*)
    integer                             :: i1, i2, i, j, m, mj
    real(kind=dbl),       allocatable   :: pmm(:), pmj0(:), pmj1(:), pmj2(:), csx(:), snx(:)
    complex(kind=dbl),    allocatable   :: ssm(:), asm(:), sumN(:), sumS(:)
    
    !Allocating memory :: step is set to 4
    allocate( pmm(4), pmj0(4), pmj1(4), pmj2(4), csx(4), snx(4), ssm(4), asm(4), &
            & sumN(4*(this%jmax+1)), sumS(4*(this%jmax+1))                       )
    
    !Cycle over latitudes
    do i = 1, this%nLegendre, 4
      do concurrent ( i2 = 1:4 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 4*(this%jmax+1), sumN(1) )
      call zero_carray_sub( 4*(this%jmax+1), sumS(1) )
      
      !Northern hemisphere :: fft towards space
      call this%fourtrans%exec_r2c_sub( 4, grid(i:i+3,:,1), sumN(1) )
      
      !Southern hemisphere :: fft towards space
      call this%fourtrans%exec_r2c_sub( 4, grid(i:i+3,:,2), sumS(1) )
      
      do m = 0, this%jmax
        call pmj_forward_rcb_4_sub( 1, this%weight(i), sumN(4*m+1), sumS(4*m+1), ssm(1), asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_4_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( 1, pmj0(1), ssm(1), cr(mj) )
        
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_4_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( 1, pmj0(1), asm(1), cr(mj-1) )
          
          call pmj_recursion_4_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( 1, pmj0(1), ssm(1), cr(mj) )
        end do
        
        if ( mod(this%jmax-m,2) /= 0 ) then
          call pmj_recursion_4_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( 1, pmj0(1), asm(1), cr(mj+1) )
        end if
      end do
    end do
    
    !Cleaning
    deallocate( pmm, pmj0, pmj1, pmj2, csx, snx, ssm, asm, sumN, sumS )
    
  end subroutine grid_to_space_sub
  
end submodule sgs