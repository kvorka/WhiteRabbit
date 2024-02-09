submodule (SphericalHarmonics) vcvv_vcvgv
  implicit none ; contains
  
  module pure subroutine vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(*), q(*), v(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    integer                           :: i, j, m, l, ijm, ijml, mj, mj1, mj2, i1, i2
    real(kind=dbl),       allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), weight(:), grid(:)
    complex(kind=dbl)                 :: cr12
    complex(kind=dbl),    allocatable :: cab(:), cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    allocate( cab(5*this%jmv1), cc(15*this%jms2), cr(4*this%jms1) )
      
      cab = czero
      cc  = czero
      cr  = czero
      
      do concurrent( ijml = 0:this%jmv-1 )
        cab(1+5*ijml) = q(ijml+1)
        cab(2+5*ijml) = v(ijml+1)
      end do
      
      call gradvec_to_vectors_sub( this%jmax, 3, 5, ri, v(1), dv_r(1), cab(1) )
      call vectors_to_scalars_sub( this%jmax, 5, cab(1), cc(1) )
      
    deallocate(cab)
    
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), weight(step), sinx(step), ssym(step*15),   &
            & asym(step*15), sumN(15*step*(this%maxj+1)), sumS(15*step*(this%maxj+1)), grid(15*step*this%nFourier) )
      
      do i = 1, this%nLegendre, step
        call lege_setup_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 15*step*(this%maxj+1), sumN(1), sumS(1) )
        
        do m = 0, this%maxj
          call pmm_recursion_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm)) < this%tolm) exit
          
          call zero_poly_sub( 15*step, ssym(1), asym(1) )
          
          j = m
            mj = m*(this%maxj+1)-m*(m+1)/2+j
            
            call pmj_setup_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 15, pmj(1), cc(1+15*mj), ssym(1) )
          
          do j = 1, (this%maxj-m)/2
            mj = mj+2
            
            call pmj_recursion_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 15, pmj(1), cc(1+15*(mj-1)), asym(1) )
            
            call pmj_recursion_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 15, pmj(1), cc(1+15*mj), ssym(1) )
            
            if ( maxval(abs(pmj)) < this%tolm ) exit
          end do
          
          if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
            mj = mj+1
            
            call pmj_recursion_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 15, pmj(1), cc(1+15*mj), asym(1) )
          end if
          
          call pmj_backward_recomb_sub( m, 15, ssym(1), asym(1), sumN(1+15*step*m), sumS(1+15*step*m) )
        end do
        
        call this%fourtrans%exec_c2r_sub(15*step, this%maxj+1, sumN, grid)
        
        do i1 = 0, this%nFourier-1
          do i2 = 0, step-1
            grid(1+4*(i2+i1*step)) = sum( grid(1+15*(i2+i1*step):3+15*(i2+i1*step)) * grid( 4+15*(i2+i1*step): 6+15*(i2+i1*step)) )
            grid(2+4*(i2+i1*step)) = sum( grid(4+15*(i2+i1*step):6+15*(i2+i1*step)) * grid( 7+15*(i2+i1*step): 9+15*(i2+i1*step)) )
            grid(3+4*(i2+i1*step)) = sum( grid(4+15*(i2+i1*step):6+15*(i2+i1*step)) * grid(10+15*(i2+i1*step):12+15*(i2+i1*step)) )
            grid(4+4*(i2+i1*step)) = sum( grid(4+15*(i2+i1*step):6+15*(i2+i1*step)) * grid(13+15*(i2+i1*step):15+15*(i2+i1*step)) )
          end do
        end do
        
        call this%fourtrans%exec_r2c_sub(4*step, this%maxj, grid, sumN)
        
        call this%fourtrans%exec_c2r_sub(15*step, this%maxj+1, sumS, grid)
        
        do i1 = 0, this%nFourier-1
          do i2 = 0, step-1
            grid(1+4*(i2+i1*step)) = sum( grid(1+15*(i2+i1*step):3+15*(i2+i1*step)) * grid( 4+15*(i2+i1*step): 6+15*(i2+i1*step)) )
            grid(2+4*(i2+i1*step)) = sum( grid(4+15*(i2+i1*step):6+15*(i2+i1*step)) * grid( 7+15*(i2+i1*step): 9+15*(i2+i1*step)) )
            grid(3+4*(i2+i1*step)) = sum( grid(4+15*(i2+i1*step):6+15*(i2+i1*step)) * grid(10+15*(i2+i1*step):12+15*(i2+i1*step)) )
            grid(4+4*(i2+i1*step)) = sum( grid(4+15*(i2+i1*step):6+15*(i2+i1*step)) * grid(13+15*(i2+i1*step):15+15*(i2+i1*step)) )
          end do
        end do
        
        call this%fourtrans%exec_r2c_sub(4*step, this%maxj, grid, sumS)
        
        do m = 0, this%jmax+1
          call pmm_recursion_sub( m, sinx(1), pmm(1) )
          call pmj_forward_recomb_sub( m, 4, weight(1), sumN(1+4*step*m), sumS(1+4*step*m), ssym(1), asym(1) )
          
          j = m
            mj = m*this%maxj-m*(m+1)/2+j
            
            call pmj_setup_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_forward_sub( 4, pmj(1), ssym(1), cr(1+4*mj) )
          
          do j = 1, (this%jmax+1-m)/2
            mj = mj+2
            
            call pmj_recursion_sub( this%amjrr(mj+m), this%bmjrr(mj+m), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_forward_sub( 4, pmj(1), asym(1), cr(1+4*(mj-1)) )
            
            call pmj_recursion_sub( this%amjrr(mj+m+1), this%bmjrr(mj+m+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_forward_sub( 4, pmj(1), ssym(1), cr(1+4*mj) )
          end do
          
          if (mod(this%jmax+1-m,2) /= 0) then
            mj = mj+1
            
            call pmj_recursion_sub( this%amjrr(mj+m+1), this%bmjrr(mj+m+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_forward_sub( 4, pmj(1), asym(1), cr(1+4*mj) )
          end if
        end do
      end do
      
    deallocate( cc, sumN, sumS, grid, pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym )
      
      do concurrent ( mj = 0:this%jms1-1 )
        cr12       = +( cr(2+4*mj) - cr(3+4*mj) * cunit ) * sq2_1
        cr(3+4*mj) = -( cr(2+4*mj) + cr(3+4*mj) * cunit ) * sq2_1
        cr(2+4*mj) = cr12
      end do
      
      j = 0
        m = 0
          ijm = 0
          mj  = m*(this%maxj)-m*(m+1)/2+j
          
          mj2 = 4*(mj + this%maxj - m)
          mj  = 4*(mj)
          
          cjm(1) = cr(1+mj )
          cjm(2) = czero
          cjm(3) = czero
          cjm(4) = cr(2+mj2) * cleb_fn(+1,-1,j,m) + cr(8+mj ) * cleb_fn(+1, 0,j,m) + conjg( cr(2+mj2) ) * cleb_fn(+1,+1,j,m)
          
          cjm(1)%im = 0._dbl
          cjm(4)%im = 0._dbl
          
      do j = 1, this%jmax
        m = 0
          ijm = ijm+4
          mj  = m*(this%maxj)-m*(m+1)/2+j
          
          mj2 = 4*(mj + this%maxj - m - 1)
          mj  = 4*(mj)
          
          cjm(1+ijm) = cr( 1+mj )
          cjm(2+ijm) = cr(-2+mj2) * cleb_fn(-1,-1,j,m) + cr(  mj) * cleb_fn(-1, 0,j,m) + conjg( cr(-2+mj2) ) * cleb_fn(-1,+1,j,m)
          cjm(3+ijm) = cr( 2+mj2) * cleb_fn( 0,-1,j,m) + cr(4+mj) * cleb_fn( 0, 0,j,m) + conjg( cr( 2+mj2) ) * cleb_fn( 0,+1,j,m)
          cjm(4+ijm) = cr( 6+mj2) * cleb_fn(+1,-1,j,m) + cr(8+mj) * cleb_fn(+1, 0,j,m) + conjg( cr( 6+mj2) ) * cleb_fn(+1,+1,j,m)
          
          cjm(1+ijm)%im = 0._dbl
          cjm(2+ijm)%im = 0._dbl
          cjm(3+ijm)%re = 0._dbl
          cjm(4+ijm)%im = 0._dbl
          
        do m = 1, j
          ijm = ijm+4
          mj  = m*(this%maxj)-m*(m+1)/2+j
          
          mj1 = 4*(mj - this%maxj + m)
          mj2 = 4*(mj + this%maxj - m - 1)
          mj  = 4*(mj)
          
          cjm(1+ijm) = cr( 1+mj )
          cjm(2+ijm) = cr(-2+mj2) * cleb_fn(-1,-1,j,m) + cr(  mj) * cleb_fn(-1, 0,j,m) + cr(-1+mj1) * cleb_fn(-1,+1,j,m)
          cjm(3+ijm) = cr( 2+mj2) * cleb_fn( 0,-1,j,m) + cr(4+mj) * cleb_fn( 0, 0,j,m) + cr( 3+mj1) * cleb_fn( 0,+1,j,m)
          cjm(4+ijm) = cr( 6+mj2) * cleb_fn(+1,-1,j,m) + cr(8+mj) * cleb_fn(+1, 0,j,m) + cr( 7+mj1) * cleb_fn(+1,+1,j,m)
        end do
      end do
      
    deallocate(cr)
    
    do concurrent ( i1 = 1:4*this%jms )
      cjm(i1) = cjm(i1) * this%scale
    end do
    
  end subroutine vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv