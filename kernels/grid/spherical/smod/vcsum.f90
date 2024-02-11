submodule (SphericalHarmonics) vcsum
  implicit none ; contains
  
  module pure subroutine vcsum_sub(this, cajm, cbjm, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    integer                           :: i, j, m, mj, i1, i2
    real(kind=dbl),       allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:), sinx(:), grid(:)
    complex(kind=dbl),    allocatable :: cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    allocate( cc(2*this%jms2), cr(this%jms1) )
      
      cc = czero
      cr = czero
      
      do m = 0, this%jmax
        do j = m, this%jmax
          cc(1+2*(m*(this%maxj+1)-m*(m+1)/2+j)) = cajm(j*(j+1)/2+m+1)
          cc(2+2*(m*(this%maxj+1)-m*(m+1)/2+j)) = cbjm(j*(j+1)/2+m+1)
        end do
      end do
    
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), weight(step), sinx(step), ssym(2*step), asym(2*step), &
            & sumN(2*step*(1+this%maxj)), sumS(2*step*(1+this%maxj)), grid(2*step*this%nFourier)                              )
      
      do i = 1, this%nLegendre, step
        call lege_setup_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 2*step*(this%maxj+1), sumN(1), sumS(1) )
        
        do m = 0, this%maxj
          call pmm_recursion_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm)) < this%tolm) exit
          
          call zero_poly_sub( 2*step, ssym(1), asym(1) )
          
          j = m
            mj = m*(this%maxj+1)-m*(m+1)/2+j
            
            call pmj_setup_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 2, pmj(1), cc(1+2*mj), ssym(1) )
          
          do j = 1, (this%maxj-m)/2
            mj = mj+2
            
            call pmj_recursion_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 2, pmj(1), cc(1+2*(mj-1)), asym(1) )
            
            call pmj_recursion_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 2, pmj(1), cc(1+2*mj), ssym(1) )
            
            if ( maxval(abs(pmj)) < this%tolm ) exit
          end do
          
          if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
            mj = mj+1
            
            call pmj_recursion_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
            call pmj_backward_sub( 2, pmj(1), cc(1+2*mj), asym(1) )
          end if
          
          call pmj_backward_recomb_sub( m, 2, ssym(1), asym(1), sumN(1+2*step*m), sumS(1+2*step*m) )
        end do
              
        call this%fourtrans%exec_c2r_sub(2*step, this%maxj+1, sumN, grid)
        
        do i1 = 0, this%nFourier-1
          do i2 = 0, step-1
            grid(1+(i2+i1*step)) = grid(1+2*(i2+i1*step)) * grid(2+2*(i2+i1*step))
          end do
        end do
        
        call this%fourtrans%exec_r2c_sub(step, this%maxj, grid, sumN)
        
        call this%fourtrans%exec_c2r_sub(2*step, this%maxj+1, sumS, grid)
        
        do i1 = 0, this%nFourier-1
          do i2 = 0, step-1
            grid(1+(i2+i1*step)) = grid(1+2*(i2+i1*step)) * grid(2+2*(i2+i1*step))
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
    
  end subroutine vcsum_sub
  
end submodule vcsum