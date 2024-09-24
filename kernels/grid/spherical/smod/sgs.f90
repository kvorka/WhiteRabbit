submodule (SphericalHarmonics) sgs
  implicit none; contains
  
  module pure subroutine space_to_grid_sub(this, cc, grid)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cc(*)
    real(kind=dbl),       intent(out) :: grid(:,:,:)
    integer                           :: i1, i2, i, j, m, mj
    real(kind=dbl),       allocatable :: pmm(:), pmj0(:), pmj1(:), pmj2(:), csx(:), snx(:)
    complex(kind=dbl),    allocatable :: ssm(:), asm(:), sumN(:), sumS(:)
    
    !Allocating memory
    allocate( pmm(16), pmj0(16), pmj1(16), pmj2(16), csx(16), snx(16),         &
            & ssm(16), asm(16), sumN(16*(this%jmax+1)), sumS(16*(this%jmax+1)) )
    
    !Cycle over latitudes :: calculating 16 at once
    do i = 1, (this%nLegendre/16)*16, 16
      do concurrent ( i2 = 1:16 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 16*(this%jmax+1), sumN(1) )
      call zero_carray_sub( 16*(this%jmax+1), sumS(1) )
      
      do m = 0, this%jmax
        call zero_carray_sub( 16, ssm(1) )
        call zero_carray_sub( 16, asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_16_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( 1, pmj0(1), cc(mj), ssm(1) )
          
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_16_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( 1, pmj0(1), cc(mj-1), asm(1) )
          
          call pmj_recursion_16_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( 1, pmj0(1), cc(mj), ssm(1) )
        end do
        
        if ( mod((this%jmax-m),2) /= 0 ) then
          call pmj_recursion_16_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( 1, pmj0(1), cc(mj+1), asm(1) )
        end if
        
        call pmj_backward_rcb_16_sub( 1, ssm(1), asm(1), sumN(16*m+1), sumS(16*m+1) )
      end do
      
      call this%fourtrans%exec_c2r_sub( 16, sumN(1), grid(i:i+15,:,1) )
      call this%fourtrans%exec_c2r_sub( 16, sumS(1), grid(i:i+15,:,2) )
    end do
    
    !Cycle over latitudes :: calculating 8 at once
    do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
      do concurrent ( i2 = 1:8 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 8*(this%jmax+1), sumN(1) )
      call zero_carray_sub( 8*(this%jmax+1), sumS(1) )
      
      do m = 0, this%jmax
        call zero_carray_sub( 8, ssm(1) )
        call zero_carray_sub( 8, asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_8_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( 1, pmj0(1), cc(mj), ssm(1) )
          
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_8_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( 1, pmj0(1), cc(mj-1), asm(1) )
          
          call pmj_recursion_8_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( 1, pmj0(1), cc(mj), ssm(1) )
        end do
        
        if ( mod((this%jmax-m),2) /= 0 ) then
          call pmj_recursion_8_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( 1, pmj0(1), cc(mj+1), asm(1) )
        end if
        
        call pmj_backward_rcb_8_sub( 1, ssm(1), asm(1), sumN(8*m+1), sumS(8*m+1) )
      end do
      
      call this%fourtrans%exec_c2r_sub( 8, sumN(1), grid(i:i+7,:,1) )
      call this%fourtrans%exec_c2r_sub( 8, sumS(1), grid(i:i+7,:,2) )
    end do
    
    !Cycle over latitudes :: calculating 4 at once
    do i = (this%nLegendre/8)*8+1, this%nLegendre, 4
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
      
      call this%fourtrans%exec_c2r_sub( 4, sumN(1), grid(i:i+3,:,1) )
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
    
    !Allocating memory
    allocate( pmm(16), pmj0(16), pmj1(16), pmj2(16), csx(16), snx(16),         &
            & ssm(16), asm(16), sumN(16*(this%jmax+1)), sumS(16*(this%jmax+1)) )
    
    !Cycle over latitudes :: computing 16 at once
    do i = 1, (this%nLegendre/16)*16, 16
      do concurrent ( i2 = 1:16 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 16*(this%jmax+1), sumN(1) )
      call zero_carray_sub( 16*(this%jmax+1), sumS(1) )
      
      call this%fourtrans%exec_r2c_sub( 16, grid(i:i+15,:,1), sumN(1) )
      call this%fourtrans%exec_r2c_sub( 16, grid(i:i+15,:,2), sumS(1) )
      
      do m = 0, this%jmax
        call pmj_forward_rcb_16_sub( 1, this%weight(i), sumN(16*m+1), sumS(16*m+1), ssm(1), asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_16_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( 1, pmj0(1), ssm(1), cr(mj) )
        
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_16_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( 1, pmj0(1), asm(1), cr(mj-1) )
          
          call pmj_recursion_16_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( 1, pmj0(1), ssm(1), cr(mj) )
        end do
        
        if ( mod(this%jmax-m,2) /= 0 ) then
          call pmj_recursion_16_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( 1, pmj0(1), asm(1), cr(mj+1) )
        end if
      end do
    end do
    
    !Cycle over latitudes :: computing 8 at once
    do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
      do concurrent ( i2 = 1:8 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 8*(this%jmax+1), sumN(1) )
      call zero_carray_sub( 8*(this%jmax+1), sumS(1) )
      
      call this%fourtrans%exec_r2c_sub( 8, grid(i:i+7,:,1), sumN(1) )
      call this%fourtrans%exec_r2c_sub( 8, grid(i:i+7,:,2), sumS(1) )
      
      do m = 0, this%jmax
        call pmj_forward_rcb_8_sub( 1, this%weight(i), sumN(8*m+1), sumS(8*m+1), ssm(1), asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_8_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( 1, pmj0(1), ssm(1), cr(mj) )
        
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_8_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( 1, pmj0(1), asm(1), cr(mj-1) )
          
          call pmj_recursion_8_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( 1, pmj0(1), ssm(1), cr(mj) )
        end do
        
        if ( mod(this%jmax-m,2) /= 0 ) then
          call pmj_recursion_8_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( 1, pmj0(1), asm(1), cr(mj+1) )
        end if
      end do
    end do
    
    !Cycle over latitudes :: computing 4 at once
    do i = (this%nLegendre/8)*8+1, this%nLegendre, 4
      do concurrent ( i2 = 1:4 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 4*(this%jmax+1), sumN(1) )
      call zero_carray_sub( 4*(this%jmax+1), sumS(1) )
      
      call this%fourtrans%exec_r2c_sub( 4, grid(i:i+3,:,1), sumN(1) )
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