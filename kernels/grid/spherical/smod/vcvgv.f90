submodule (SphericalHarmonics) vcvgv
  implicit none ; contains
  
  module pure subroutine vcvgv_sub(this, ri, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(*), v(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    integer                           :: i, j, m, ijm, mj, mj1, mj2, lmj, i1, i2
    real(kind=dbl),       allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:), sinx(:), grid(:,:,:), fft(:,:,:)
    complex(kind=dbl),    allocatable :: cab(:), cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    allocate( cab(4*this%jmv1), cc(12*this%jms2), cr(3*this%jms1) )
      
      cab = czero
      cc  = czero
      cr  = czero
      
      do concurrent( ijm = 0:this%jmv-1 )
        cab(1+4*ijm) = v(ijm+1)
      end do
      
      call gradvec_to_vectors_sub( this%jmax, 2, 4, ri, v(1), dv_r(1), cab(1) )
      call vectors_to_scalars_sub( this%jmax, 4, cab(1), cc(1) )
      
    deallocate(cab)
    
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), weight(step), sinx(step), ssym(step*12), asym(step*12), &
            & sumN(12*step*(1+this%maxj)), sumS(12*step*(1+this%maxj)), grid(12,step,0:this%nFourier-1), &
            & fft(3,step,0:this%nFourier-1) )
    
    do i = 1, this%nLegendre, step
      call lege_setup_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
      call zero_poly_sub( 12*step*(this%maxj+1), sumN(1), sumS(1) )
      
      do m = 0, this%maxj
        call pmm_recursion_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm)) < this%tolm) exit
        
        call zero_poly_sub( 12*step, ssym(1), asym(1) )
        
        j = m
          mj = m*(this%maxj+1)-m*(m+1)/2+j
          
          call pmj_setup_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
          call pmj_backward_sub( 12, pmj(1), cc(1+12*mj), ssym(1) )
        
        do j = 1, (this%maxj-m)/2
          mj = mj+2
          
          call pmj_recursion_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
          call pmj_backward_sub( 12, pmj(1), cc(1+12*(mj-1)), asym(1) )
          
          call pmj_recursion_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
          call pmj_backward_sub( 12, pmj(1), cc(1+12*mj), ssym(1) )
          
          if ( maxval(abs(pmj)) < this%tolm ) exit
        end do
        
        if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
          mj = mj+1
          
          call pmj_recursion_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
          call pmj_backward_sub( 12, pmj(1), cc(1+12*mj), asym(1) )
        end if
        
        call pmj_backward_recomb_sub( m, 12, ssym(1), asym(1), sumN(1+12*step*m), sumS(1+12*step*m) )
      end do
      
      call this%fourtrans%exec_c2r_sub(12*step, this%maxj+1, sumN, grid)
      
      do concurrent ( i1=0:this%nFourier-1, i2=1:step )
        fft(1,i2,i1) = sum( grid(1:3,i2,i1) * grid( 4: 6,i2,i1) )
        fft(2,i2,i1) = sum( grid(1:3,i2,i1) * grid( 7: 9,i2,i1) )
        fft(3,i2,i1) = sum( grid(1:3,i2,i1) * grid(10:12,i2,i1) )
      end do
      
      call this%fourtrans%exec_r2c_sub(3*step, this%maxj, fft, sumN)
      
      call this%fourtrans%exec_c2r_sub(12*step, this%maxj+1, sumS, grid)
      
      do concurrent ( i1=0:this%nFourier-1, i2=1:step )
        fft(1,i2,i1) = sum( grid(1:3,i2,i1) * grid( 4: 6,i2,i1) )
        fft(2,i2,i1) = sum( grid(1:3,i2,i1) * grid( 7: 9,i2,i1) )
        fft(3,i2,i1) = sum( grid(1:3,i2,i1) * grid(10:12,i2,i1) )
      end do
      
      call this%fourtrans%exec_r2c_sub(3*step, this%maxj, fft, sumS)
      
      do m = 0, this%jmax+1
        call pmm_recursion_sub( m, sinx(1), pmm(1) )
        call pmj_forward_recomb_sub( m, 3, weight(1), sumN(1+3*step*m), sumS(1+3*step*m), ssym(1), asym(1) )
        
        j = m
          mj = m*this%maxj-m*(m+1)/2+j
          
          call pmj_setup_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
          call pmj_forward_sub( 3, pmj(1), ssym(1), cr(1+3*mj) )
        
        do j = 1, (this%jmax+1-m)/2
          mj = mj+2
          
          call pmj_recursion_sub( this%amjrr(mj+m), this%bmjrr(mj+m), cosx(1), pmj2(1), pmj1(1), pmj(1) )
          call pmj_forward_sub( 3, pmj(1), asym(1), cr(1+3*(mj-1)) )
          
          call pmj_recursion_sub( this%amjrr(mj+m+1), this%bmjrr(mj+m+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
          call pmj_forward_sub( 3, pmj(1), ssym(1), cr(1+3*mj) )
        end do
        
        if (mod(this%jmax+1-m,2) /= 0) then
          mj = mj+1
          
          call pmj_recursion_sub( this%amjrr(mj+m+1), this%bmjrr(mj+m+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
          call pmj_forward_sub( 3, pmj(1), asym(1), cr(1+3*mj) )
        end if
      end do
    end do
    
    deallocate( cc, sumN, sumS, grid, fft, pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym )
      
      call cartesian_to_cyclic_sub( 1, 3, this%jms1, cr(1) )
      
      j = 0
        m = 0
          ijm = 0 ; mj = m*(this%maxj)-m*(m+1)/2+j
          
          mj2 = 3*(mj + this%maxj - m)
          mj  = 3*(mj)
          
          cjm(1) = czero
          cjm(2) = czero
          cjm(3) = cr(1+mj2) * cleb_fn(+1,-1,j,m) + cr(6+mj) * cleb_fn(+1,0,j,m) + conjg( cr(1+mj2) ) * cleb_fn(+1,+1,j,m)
          
          cjm(3)%im = 0._dbl
          
      do j = 1, this%jmax
        m = 0
          ijm = ijm+3 ; mj = m*(this%maxj)-m*(m+1)/2+j
          
          mj2 = 3*(mj + this%maxj - m - 1)
          mj  = 3*(mj)
          
          cjm(1+ijm) = cr(-2+mj2) * cleb_fn(-1,-1,j,m) + cr(  mj) * cleb_fn(-1,0,j,m) + conjg( cr(-2+mj2) ) * cleb_fn(-1,+1,j,m)
          cjm(2+ijm) = cr( 1+mj2) * cleb_fn( 0,-1,j,m) + cr(3+mj) * cleb_fn( 0,0,j,m) + conjg( cr( 1+mj2) ) * cleb_fn( 0,+1,j,m)
          cjm(3+ijm) = cr( 4+mj2) * cleb_fn(+1,-1,j,m) + cr(6+mj) * cleb_fn(+1,0,j,m) + conjg( cr( 4+mj2) ) * cleb_fn(+1,+1,j,m)
          
          cjm(1+ijm)%im = 0._dbl
          cjm(2+ijm)%re = 0._dbl
          cjm(3+ijm)%im = 0._dbl
          
        do m = 1, j
          ijm = ijm+3 ; mj = m*(this%maxj)-m*(m+1)/2+j
          
          mj1 = 3*(mj - this%maxj + m)
          mj2 = 3*(mj + this%maxj - m - 1)
          mj  = 3*(mj)
          
          cjm(1+ijm) = cr(-2+mj2) * cleb_fn(-1,-1,j,m) + cr(  mj) * cleb_fn(-1,0,j,m) + cr(-1+mj1) * cleb_fn(-1,+1,j,m)
          cjm(2+ijm) = cr( 1+mj2) * cleb_fn( 0,-1,j,m) + cr(3+mj) * cleb_fn( 0,0,j,m) + cr( 2+mj1) * cleb_fn( 0,+1,j,m)
          cjm(3+ijm) = cr( 4+mj2) * cleb_fn(+1,-1,j,m) + cr(6+mj) * cleb_fn(+1,0,j,m) + cr( 5+mj1) * cleb_fn(+1,+1,j,m)
        end do
      end do
      
    deallocate(cr)
    
    do concurrent ( i1=1:3*this%jms )
      cjm(i1) = cjm(i1) * this%scale
    end do
    
  end subroutine vcvgv_sub
  
end submodule vcvgv
