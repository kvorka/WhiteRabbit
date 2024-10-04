submodule (lateral_grid) transform
  implicit none; contains
  
  module pure subroutine transform_sub(this, nf, nb, cc, cr, grid_sub)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: nf, nb
    complex(kind=dbl),    intent(in)    :: cc(nb,*)
    complex(kind=dbl),    intent(inout) :: cr(nf,*)
    integer                             :: i1, i2, i, j, m, mj
    real(kind=dbl),       allocatable   :: pmm(:), pmj0(:), pmj1(:), pmj2(:), csx(:), snx(:), grid(:)
    complex(kind=dbl),    allocatable   :: ssm(:), asm(:), sumN(:), sumS(:)
    
    interface
      module pure subroutine grid_sub(step, nfour, gxyz)
        integer,                intent(in)    :: step, nfour
        real(kind=dbl), target, intent(inout) :: gxyz(*)
      end subroutine grid_sub
    end interface
    
    !Allocating memory :: step is set to 16
    allocate( pmm(16), pmj0(16), pmj1(16), pmj2(16), csx(16), snx(16), ssm(16*nb), asm(16*nb), &
            & sumN(16*nb*(this%jmax+1)), sumS(16*nb*(this%jmax+1)), grid(16*nb*this%nFourier)  )
    
    !Cycle over latitudes :: calculating 16 at once
    do i = 1, (this%nLegendre/16)*16, 16
      do concurrent ( i2 = 1:16 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 16*nb*(this%jmax+1), sumN(1) )
      call zero_carray_sub( 16*nb*(this%jmax+1), sumS(1) )
      
      do m = 0, this%jmax
        call zero_carray_sub( 16*nb, ssm(1) )
        call zero_carray_sub( 16*nb, asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_16_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
          
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_16_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( nb, pmj0(1), cc(1,mj-1), asm(1) )
          
          call pmj_recursion_16_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
        end do
        
        if ( mod((this%jmax-m),2) /= 0 ) then
          call pmj_recursion_16_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( nb, pmj0(1), cc(1,mj+1), asm(1) )
        end if
        
        call pmj_backward_rcb_16_sub( nb, ssm(1), asm(1), sumN(16*nb*m+1), sumS(16*nb*m+1) )
      end do
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 16*nb, sumN(1), grid(1) )
      call grid_sub( 16, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 16*nf, grid(1), sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 16*nb, sumS(1), grid(1) )
      call grid_sub( 16, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 16*nf, grid(1), sumS(1) )
      
      do m = 0, this%jmax
        call pmj_forward_rcb_16_sub( nf, this%weight(i), sumN(16*nf*m+1), sumS(16*nf*m+1), ssm(1), asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_16_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_16_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( nf, pmj0(1), asm(1), cr(1,mj-1) )
          
          call pmj_recursion_16_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        end do
        
        if ( mod(this%jmax-m,2) /= 0 ) then
          call pmj_recursion_16_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( nf, pmj0(1), asm(1), cr(1,mj+1) )
        end if
      end do
    end do

    !Cycle over latitudes :: calculating 8 at once
    do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
      do concurrent ( i2 = 1:8 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 8*nb*(this%jmax+1), sumN(1) )
      call zero_carray_sub( 8*nb*(this%jmax+1), sumS(1) )
      
      do m = 0, this%jmax
        call zero_carray_sub( 8*nb, ssm(1) )
        call zero_carray_sub( 8*nb, asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_8_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
          
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_8_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( nb, pmj0(1), cc(1,mj-1), asm(1) )
          
          call pmj_recursion_8_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
        end do
        
        if ( mod((this%jmax-m),2) /= 0 ) then
          call pmj_recursion_8_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( nb, pmj0(1), cc(1,mj+1), asm(1) )
        end if
        
        call pmj_backward_rcb_8_sub( nb, ssm(1), asm(1), sumN(8*nb*m+1), sumS(8*nb*m+1) )
      end do
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 8*nb, sumN(1), grid(1) )
      call grid_sub( 8, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 8*nf, grid(1), sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 8*nb, sumS(1), grid(1) )
      call grid_sub( 8, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 8*nf, grid(1), sumS(1) )
      
      do m = 0, this%jmax
        call pmj_forward_rcb_8_sub( nf, this%weight(i), sumN(8*nf*m+1), sumS(8*nf*m+1), ssm(1), asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_8_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_8_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( nf, pmj0(1), asm(1), cr(1,mj-1) )
          
          call pmj_recursion_8_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        end do
        
        if ( mod(this%jmax-m,2) /= 0 ) then
          call pmj_recursion_8_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( nf, pmj0(1), asm(1), cr(1,mj+1) )
        end if
      end do
    end do
    
    !Cycle over latitudes :: calculating 4 possibly left at once
    do i = (this%nLegendre/8)*8+1, this%nLegendre, 4
      do concurrent ( i2 = 1:4 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 4*nb*(this%jmax+1), sumN(1) )
      call zero_carray_sub( 4*nb*(this%jmax+1), sumS(1) )
      
      do m = 0, this%jmax
        call zero_carray_sub( 4*nb, ssm(1) )
        call zero_carray_sub( 4*nb, asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_4_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
          
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_4_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( nb, pmj0(1), cc(1,mj-1), asm(1) )
          
          call pmj_recursion_4_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
        end do
        
        if ( mod((this%jmax-m),2) /= 0 ) then
          call pmj_recursion_4_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( nb, pmj0(1), cc(1,mj+1), asm(1) )
        end if
        
        call pmj_backward_rcb_4_sub( nb, ssm(1), asm(1), sumN(4*nb*m+1), sumS(4*nb*m+1) )
      end do
      
      !Northern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 4*nb, sumN(1), grid(1) )
      call grid_sub( 4, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 4*nf, grid(1), sumN(1) )
      
      !Southern hemisphere :: fft towards grid, multiplication on the grid, fft towards space
      call this%fourtrans%exec_c2r_sub( 4*nb, sumS(1), grid(1) )
      call grid_sub( 4, this%nFourier, grid(1) )
      call this%fourtrans%exec_r2c_sub( 4*nf, grid(1), sumS(1) )
      
      do m = 0, this%jmax
        call pmj_forward_rcb_4_sub( nf, this%weight(i), sumN(4*nf*m+1), sumS(4*nf*m+1), ssm(1), asm(1) )
        
        !j = m
          mj = m*(this%jmax+1)-(m-2)*(m+1)/2
          
          call pmj_mmset_4_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        
        do j = 1, (this%jmax-m)/2
          mj = mj+2
          
          call pmj_recursion_4_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( nf, pmj0(1), asm(1), cr(1,mj-1) )
          
          call pmj_recursion_4_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        end do
        
        if ( mod(this%jmax-m,2) /= 0 ) then
          call pmj_recursion_4_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( nf, pmj0(1), asm(1), cr(1,mj+1) )
        end if
      end do
    end do
    
    !Cleaning
    deallocate( pmm, pmj0, pmj1, pmj2, csx, snx, ssm, asm, sumN, sumS, grid )
    
  end subroutine transform_sub
  
end submodule transform