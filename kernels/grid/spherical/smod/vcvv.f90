submodule (SphericalHarmonics) vcvv
  implicit none ; contains
  
  module pure subroutine vcvv_sub(this, cajml, cbjml, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajml(*), cbjml(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    integer                           :: i, j, m, mj, i1, i2
    real(kind=dbl),       allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:), sinx(:), grid(:)
    complex(kind=dbl),    allocatable :: cab(:), cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    allocate( cab(2*this%jmv1), cc(6*this%jms2), cr(this%jms1) )
      
      cab = czero
      cc  = czero
      cr  = czero
      
      do concurrent( mj = 0:this%jmv-1 )
        cab(1+2*mj) = cajml(mj+1)
        cab(2+2*mj) = cbjml(mj+1)
      end do
      
      call vectors_to_scalars_sub( this%jmax, 2, cab(1), cc(1) )
      
    deallocate(cab)
    
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), weight(step), sinx(step), ssym(6*step), asym(6*step), &
            & sumN(6*step*(1+this%maxj)), sumS(6*step*(1+this%maxj)), grid(6*step*this%nFourier)                              )
      
      do i = 1, this%nLegendre, step
        call lege_setup_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 6*step*(this%maxj+1), sumN(1), sumS(1) )
        
        do m = 0, this%maxj
          call pmm_recursion_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm)) < this%tolm) exit
          
          call zero_poly_sub( 6*step, ssym(1), asym(1) )
          
          j = m
            mj = m*(this%maxj+1)-m*(m+1)/2+j
            
            call pmj_setup_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 6, pmj(1), cc(1+6*mj), ssym(1) )
          
          do j = 1, (this%maxj-m)/2
            mj = mj+2
            
            call pmj_recursion_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 6, pmj(1), cc(1+6*(mj-1)), asym(1) )
            
            call pmj_recursion_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 6, pmj(1), cc(1+6*mj), ssym(1) )
            
            if ( maxval(abs(pmj)) < this%tolm ) exit
          end do
          
          if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
            mj = mj+1
            
            call pmj_recursion_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 6, pmj(1), cc(1+6*mj), asym(1) )
          end if
          
          call pmj_backward_recomb_sub( m, 6, ssym(1), asym(1), sumN(1+6*step*m), sumS(1+6*step*m) )
        end do
              
        call this%fourtrans%exec_c2r_sub(6*step, this%maxj+1, sumN, grid)
        
        do i1 = 0, this%nFourier-1
          do i2 = 0, step-1
            grid(1+(i2+i1*step)) = sum( grid(1+6*(i2+i1*step):3+6*(i2+i1*step)) * grid( 4+6*(i2+i1*step): 6+6*(i2+i1*step)) )
          end do
        end do
        
        call this%fourtrans%exec_r2c_sub(step, this%maxj, grid, sumN)
        
        call this%fourtrans%exec_c2r_sub(6*step, this%maxj+1, sumS, grid)
        
        do i1 = 0, this%nFourier-1
          do i2 = 0, step-1
            grid(1+(i2+i1*step)) = sum( grid(1+6*(i2+i1*step):3+6*(i2+i1*step)) * grid( 4+6*(i2+i1*step): 6+6*(i2+i1*step)) )
          end do
        end do
        
        call this%fourtrans%exec_r2c_sub(step, this%maxj, grid, sumS)
        
        do m = 0, this%jmax+1
          call pmm_recursion_sub( m, sinx(1), pmm(1) )
          call pmj_forward_recomb_sub( m, 1, weight(1), sumN(1+step*m), sumS(1+step*m), ssym(1), asym(1) )
          
          j = m
            mj = m*this%maxj-m*(m+1)/2+j
            
            call pmj_setup_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_forward_sub( 1, pmj(1), ssym(1), cr(1+mj) )
          
          do j = 1, (this%jmax+1-m)/2
            mj = mj+2
            
            call pmj_recursion_sub( this%amjrr(mj+m), this%bmjrr(mj+m), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_forward_sub( 1, pmj(1), asym(1), cr(mj) )
            
            call pmj_recursion_sub( this%amjrr(mj+m+1), this%bmjrr(mj+m+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_forward_sub( 1, pmj(1), ssym(1), cr(1+mj) )
          end do
          
          if (mod(this%jmax+1-m,2) /= 0) then
            mj = mj+1
            
            call pmj_recursion_sub( this%amjrr(mj+m+1), this%bmjrr(mj+m+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_forward_sub( 1, pmj(1), asym(1), cr(1+mj) )
          end if
        end do
      end do
      
    deallocate( cc, sumN, sumS, grid, pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym )
      
      do j = 0, this%jmax
        m = 0
          cjm(j*(j+1)/2+m+1)%re = cr(m*(this%maxj)-m*(m+1)/2+j+1)%re * this%scale
          cjm(j*(j+1)/2+m+1)%im = 0._dbl
        
        do m = 1, j
          cjm(j*(j+1)/2+m+1) = cr(m*(this%maxj)-m*(m+1)/2+j+1) * this%scale
        end do
      end do
      
    deallocate(cr)
    
  end subroutine vcvv_sub
  
end submodule vcvv
