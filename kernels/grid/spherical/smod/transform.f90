submodule (SphericalHarmonics) transform
  implicit none; contains
  
  pure subroutine pmj_mmset_4_sub(m, cmm, snx, pmm, pmj2, pmj1, pmj0)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: cmm, snx(4)
    real(kind=dbl), intent(inout) :: pmm(4)
    real(kind=dbl), intent(out)   :: pmj2(4), pmj1(4), pmj0(4)
    integer                       :: i2
    
    if ( m /= 0 ) then
      do concurrent ( i2 = 1:4 )
        pmm(i2) = cmm * snx(i2) * pmm(i2)
      end do
    else
      do concurrent ( i2 = 1:4 )
        pmm(i2) = cmm
      end do
    end if
    
    do concurrent ( i2 = 1:4 )
      pmj2(i2) = zero
      pmj1(i2) = zero
      pmj0(i2) = pmm(i2)
    end do
    
  end subroutine pmj_mmset_4_sub
  
  pure subroutine pmj_recursion_4_sub(amj, bmj, csx, pmj2, pmj1, pmj0)
    real(kind=dbl), intent(in)    :: amj, bmj, csx(4)
    real(kind=dbl), intent(inout) :: pmj2(4), pmj1(4), pmj0(4)
    integer                       :: i2
    
    do concurrent ( i2=1:4 )
      pmj2(i2) = pmj1(i2)
      pmj1(i2) = pmj0(i2)
      pmj0(i2) = amj * csx(i2) * pmj1(i2) - bmj * pmj2(i2)
    end do
    
  end subroutine pmj_recursion_4_sub
  
  pure subroutine pmj_backward_sum_4_sub(nb, legep, cc, legesum)
    integer,           intent(in)    :: nb
    real(kind=dbl),    intent(in)    :: legep(4)
    complex(kind=dbl), intent(in)    :: cc(nb)
    complex(kind=dbl), intent(inout) :: legesum(4,nb)
    integer                          :: i1, i2
    
    do concurrent ( i1 = 1:nb, i2 = 1:4 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * legep(i2)
    end do
    
  end subroutine pmj_backward_sum_4_sub
  
  pure subroutine pmj_backward_rcb_4_sub(nb, sumsym, sumasym, sumN, sumS)
    integer,           intent(in)    :: nb
    complex(kind=dbl), intent(in)    :: sumsym(4,nb), sumasym(4,nb)
    complex(kind=dbl), intent(inout) :: sumN(nb,4), sumS(nb,4)
    integer                          :: i2, i1
    
    do concurrent ( i2=1:4, i1=1:nb )
      sumN(i1,i2) = sumsym(i2,i1) + sumasym(i2,i1)
      sumS(i1,i2) = sumsym(i2,i1) - sumasym(i2,i1)
    end do
    
  end subroutine pmj_backward_rcb_4_sub
  
  pure subroutine pmj_forward_sum_4_sub(nf, legep, legesum, cr)
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: legep(4)
    complex(kind=dbl), intent(in)    :: legesum(4,nf)
    complex(kind=dbl), intent(inout) :: cr(nf)
    integer                          :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + sum( legep(1:4) * legesum(1:4,i1) )
    end do
    
  end subroutine pmj_forward_sum_4_sub
  
  pure subroutine pmj_forward_rcb_4_sub(nf, w, sumN, sumS, sumsym, sumasym)
    integer,           intent(in)  :: nf
    real(kind=dbl),    intent(in)  :: w(4)
    complex(kind=dbl), intent(in)  :: sumN(nf,4), sumS(nf,4)
    complex(kind=dbl), intent(out) :: sumsym(4,nf), sumasym(4,nf)
    integer                        :: i1, i2
    
    do concurrent ( i1 = 1:nf, i2 = 1:4 )
      sumsym(i2,i1)  = w(i2) * ( sumN(i1,i2) + sumS(i1,i2) )
      sumasym(i2,i1) = w(i2) * ( sumN(i1,i2) - sumS(i1,i2) )
    end do
    
  end subroutine pmj_forward_rcb_4_sub
  
  pure subroutine pmj_mmset_8_sub(m, cmm, snx, pmm, pmj2, pmj1, pmj0)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: cmm, snx(8)
    real(kind=dbl), intent(inout) :: pmm(8)
    real(kind=dbl), intent(out)   :: pmj2(8), pmj1(8), pmj0(8)
    integer                       :: i2
    
    if ( m /= 0 ) then
      do concurrent ( i2 = 1:8 )
        pmm(i2) = cmm * snx(i2) * pmm(i2)
      end do
    else
      do concurrent ( i2 = 1:8 )
        pmm(i2) = cmm
      end do
    end if
    
    do concurrent ( i2 = 1:8 )
      pmj2(i2) = zero
      pmj1(i2) = zero
      pmj0(i2) = pmm(i2)
    end do
    
  end subroutine pmj_mmset_8_sub
  
  pure subroutine pmj_recursion_8_sub(amj, bmj, csx, pmj2, pmj1, pmj0)
    real(kind=dbl), intent(in)    :: amj, bmj, csx(8)
    real(kind=dbl), intent(inout) :: pmj2(8), pmj1(8), pmj0(8)
    integer                       :: i2
    
    do concurrent ( i2=1:8 )
      pmj2(i2) = pmj1(i2)
      pmj1(i2) = pmj0(i2)
      pmj0(i2) = amj * csx(i2) * pmj1(i2) - bmj * pmj2(i2)
    end do
    
  end subroutine pmj_recursion_8_sub
  
  pure subroutine pmj_backward_sum_8_sub(nb, legep, cc, legesum)
    integer,           intent(in)    :: nb
    real(kind=dbl),    intent(in)    :: legep(8)
    complex(kind=dbl), intent(in)    :: cc(nb)
    complex(kind=dbl), intent(inout) :: legesum(8,nb)
    integer                          :: i1, i2
    
    do concurrent ( i1 = 1:nb, i2 = 1:8 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * legep(i2)
    end do
    
  end subroutine pmj_backward_sum_8_sub
  
  pure subroutine pmj_backward_rcb_8_sub(nb, sumsym, sumasym, sumN, sumS)
    integer,           intent(in)    :: nb
    complex(kind=dbl), intent(in)    :: sumsym(8,nb), sumasym(8,nb)
    complex(kind=dbl), intent(inout) :: sumN(nb,8), sumS(nb,8)
    integer                          :: i2, i1
    
    do concurrent ( i2=1:8, i1=1:nb )
      sumN(i1,i2) = sumsym(i2,i1) + sumasym(i2,i1)
      sumS(i1,i2) = sumsym(i2,i1) - sumasym(i2,i1)
    end do
    
  end subroutine pmj_backward_rcb_8_sub
  
  pure subroutine pmj_forward_sum_8_sub(nf, legep, legesum, cr)
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: legep(8)
    complex(kind=dbl), intent(in)    :: legesum(8,nf)
    complex(kind=dbl), intent(inout) :: cr(nf)
    integer                          :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + sum( legep(1:8) * legesum(1:8,i1) )
    end do
    
  end subroutine pmj_forward_sum_8_sub
  
  pure subroutine pmj_forward_rcb_8_sub(nf, w, sumN, sumS, sumsym, sumasym)
    integer,           intent(in)  :: nf
    real(kind=dbl),    intent(in)  :: w(8)
    complex(kind=dbl), intent(in)  :: sumN(nf,8), sumS(nf,8)
    complex(kind=dbl), intent(out) :: sumsym(8,nf), sumasym(8,nf)
    integer                        :: i1, i2
    
    do concurrent ( i1 = 1:nf, i2 = 1:8 )
      sumsym(i2,i1)  = w(i2) * ( sumN(i1,i2) + sumS(i1,i2) )
      sumasym(i2,i1) = w(i2) * ( sumN(i1,i2) - sumS(i1,i2) )
    end do
    
  end subroutine pmj_forward_rcb_8_sub
  
  pure subroutine pmj_mmset_16_sub(m, cmm, snx, pmm, pmj2, pmj1, pmj0)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: cmm, snx(16)
    real(kind=dbl), intent(inout) :: pmm(16)
    real(kind=dbl), intent(out)   :: pmj2(16), pmj1(16), pmj0(16)
    integer                       :: i2
    
    if ( m /= 0 ) then
      do concurrent ( i2 = 1:16 )
        pmm(i2) = cmm * snx(i2) * pmm(i2)
      end do
    else
      do concurrent ( i2 = 1:16 )
        pmm(i2) = cmm
      end do
    end if
    
    do concurrent ( i2 = 1:16 )
      pmj2(i2) = zero
      pmj1(i2) = zero
      pmj0(i2) = pmm(i2)
    end do
    
  end subroutine pmj_mmset_16_sub
  
  pure subroutine pmj_recursion_16_sub(amj, bmj, csx, pmj2, pmj1, pmj0)
    real(kind=dbl), intent(in)    :: amj, bmj, csx(16)
    real(kind=dbl), intent(inout) :: pmj2(16), pmj1(16), pmj0(16)
    integer                       :: i2
    
    do concurrent ( i2=1:16 )
      pmj2(i2) = pmj1(i2)
      pmj1(i2) = pmj0(i2)
      pmj0(i2) = amj * csx(i2) * pmj1(i2) - bmj * pmj2(i2)
    end do
    
  end subroutine pmj_recursion_16_sub
  
  pure subroutine pmj_backward_sum_16_sub(nb, legep, cc, legesum)
    integer,           intent(in)    :: nb
    real(kind=dbl),    intent(in)    :: legep(16)
    complex(kind=dbl), intent(in)    :: cc(nb)
    complex(kind=dbl), intent(inout) :: legesum(16,nb)
    integer                          :: i1, i2
    
    do concurrent ( i1 = 1:nb, i2 = 1:16 )
      legesum(i2,i1) = legesum(i2,i1) + cc(i1) * legep(i2)
    end do
    
  end subroutine pmj_backward_sum_16_sub
  
  pure subroutine pmj_backward_rcb_16_sub(nb, sumsym, sumasym, sumN, sumS)
    integer,           intent(in)    :: nb
    complex(kind=dbl), intent(in)    :: sumsym(16,nb), sumasym(16,nb)
    complex(kind=dbl), intent(inout) :: sumN(nb,16), sumS(nb,16)
    integer                          :: i2, i1
    
    do concurrent ( i2=1:16, i1=1:nb )
      sumN(i1,i2) = sumsym(i2,i1) + sumasym(i2,i1)
      sumS(i1,i2) = sumsym(i2,i1) - sumasym(i2,i1)
    end do
    
  end subroutine pmj_backward_rcb_16_sub
  
  pure subroutine pmj_forward_sum_16_sub(nf, legep, legesum, cr)
    integer,           intent(in)    :: nf
    real(kind=dbl),    intent(in)    :: legep(16)
    complex(kind=dbl), intent(in)    :: legesum(16,nf)
    complex(kind=dbl), intent(inout) :: cr(nf)
    integer                          :: i1
    
    do concurrent ( i1=1:nf )
      cr(i1) = cr(i1) + sum( legep(1:16) * legesum(1:16,i1) )
    end do
    
  end subroutine pmj_forward_sum_16_sub
  
  pure subroutine pmj_forward_rcb_16_sub(nf, w, sumN, sumS, sumsym, sumasym)
    integer,           intent(in)  :: nf
    real(kind=dbl),    intent(in)  :: w(16)
    complex(kind=dbl), intent(in)  :: sumN(nf,16), sumS(nf,16)
    complex(kind=dbl), intent(out) :: sumsym(16,nf), sumasym(16,nf)
    integer                        :: i1, i2
    
    do concurrent ( i1 = 1:nf, i2 = 1:16 )
      sumsym(i2,i1)  = w(i2) * ( sumN(i1,i2) + sumS(i1,i2) )
      sumasym(i2,i1) = w(i2) * ( sumN(i1,i2) - sumS(i1,i2) )
    end do
    
  end subroutine pmj_forward_rcb_16_sub
  
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
            & sumN(16*nb*this%jmax3), sumS(16*nb*this%jmax3), grid(16*nb*this%nFourier)        )
    
    !Cycle over latitudes :: calculating 16 at once
    do i = 1, (this%nLegendre/16)*16, 16
      do concurrent ( i2 = 1:16 )
        csx(i2) = this%cosx(i+i2-1)
        snx(i2) = sqrt(1-csx(i2)**2)
      end do
      
      call zero_carray_sub( 16*nb*this%jmax3, sumN(1) )
      call zero_carray_sub( 16*nb*this%jmax3, sumS(1) )
      
      do m = 0, this%jmax2
        call zero_carray_sub( 16*nb, ssm(1) )
        call zero_carray_sub( 16*nb, asm(1) )
        
        !j = m
          mj = m*this%jmax3-(m-2)*(m+1)/2
          
          call pmj_mmset_16_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call pmj_recursion_16_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( nb, pmj0(1), cc(1,mj-1), asm(1) )
          
          call pmj_recursion_16_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_16_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
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
      
      do m = 0, this%jmax2
        call pmj_forward_rcb_16_sub( nf, this%weight(i), sumN(16*nf*m+1), sumS(16*nf*m+1), ssm(1), asm(1) )
        
        !j = m
          mj = m*this%jmax3-(m-2)*(m+1)/2
          
          call pmj_mmset_16_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call pmj_recursion_16_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( nf, pmj0(1), asm(1), cr(1,mj-1) )
          
          call pmj_recursion_16_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_16_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
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
      
      call zero_carray_sub( 8*nb*this%jmax3, sumN(1) )
      call zero_carray_sub( 8*nb*this%jmax3, sumS(1) )
      
      do m = 0, this%jmax2
        call zero_carray_sub( 8*nb, ssm(1) )
        call zero_carray_sub( 8*nb, asm(1) )
        
        !j = m
          mj = m*this%jmax3-(m-2)*(m+1)/2
          
          call pmj_mmset_8_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call pmj_recursion_8_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( nb, pmj0(1), cc(1,mj-1), asm(1) )
          
          call pmj_recursion_8_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_8_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
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
      
      do m = 0, this%jmax2
        call pmj_forward_rcb_8_sub( nf, this%weight(i), sumN(8*nf*m+1), sumS(8*nf*m+1), ssm(1), asm(1) )
        
        !j = m
          mj = m*this%jmax3-(m-2)*(m+1)/2
          
          call pmj_mmset_8_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call pmj_recursion_8_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( nf, pmj0(1), asm(1), cr(1,mj-1) )
          
          call pmj_recursion_8_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_8_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
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
      
      call zero_carray_sub( 4*nb*this%jmax3, sumN(1) )
      call zero_carray_sub( 4*nb*this%jmax3, sumS(1) )
      
      do m = 0, this%jmax2
        call zero_carray_sub( 4*nb, ssm(1) )
        call zero_carray_sub( 4*nb, asm(1) )
        
        !j = m
          mj = m*this%jmax3-(m-2)*(m+1)/2
          
          call pmj_mmset_4_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
          
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call pmj_recursion_4_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( nb, pmj0(1), cc(1,mj-1), asm(1) )
          
          call pmj_recursion_4_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_backward_sum_4_sub( nb, pmj0(1), cc(1,mj), ssm(1) )
        end do
        
        if ( mod((this%jmax2-m),2) /= 0 ) then
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
      
      do m = 0, this%jmax2
        call pmj_forward_rcb_4_sub( nf, this%weight(i), sumN(4*nf*m+1), sumS(4*nf*m+1), ssm(1), asm(1) )
        
        !j = m
          mj = m*this%jmax3-(m-2)*(m+1)/2
          
          call pmj_mmset_4_sub( m, this%cmm(m), snx(1), pmm(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        
        do j = 1, (this%jmax2-m)/2
          mj = mj+2
          
          call pmj_recursion_4_sub( this%amj(mj-1), this%bmj(mj-1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( nf, pmj0(1), asm(1), cr(1,mj-1) )
          
          call pmj_recursion_4_sub( this%amj(mj), this%bmj(mj), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( nf, pmj0(1), ssm(1), cr(1,mj) )
        end do
        
        if ( mod(this%jmax2-m,2) /= 0 ) then
          call pmj_recursion_4_sub( this%amj(mj+1), this%bmj(mj+1), csx(1), pmj2(1), pmj1(1), pmj0(1) )
          call pmj_forward_sum_4_sub( nf, pmj0(1), asm(1), cr(1,mj+1) )
        end if
      end do
    end do
    
    !Cleaning
    deallocate( pmm, pmj0, pmj1, pmj2, csx, snx, ssm, asm, sumN, sumS, grid )
    
  end subroutine transform_sub
  
end submodule transform