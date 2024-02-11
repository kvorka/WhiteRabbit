submodule (SphericalHarmonics) vctol
  implicit none ; contains
  
  module subroutine vctol_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer                             :: i, j, m, mj, i1, i2
    real(kind=dbl)                      :: diff, maxdiff
    real(kind=dbl),       allocatable   :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:), sinx(:)
    complex(kind=dbl),    allocatable   :: cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    allocate( cc(this%jms2), cr(this%jms1) )
      
      cc = czero
      
      do m = 0, this%jmax
        do j = m, this%jmax
          if ( m /= 0 ) then
            call random_number(cc(m*(this%maxj+1)-m*(m+1)/2+j+1)%re)
            call random_number(cc(m*(this%maxj+1)-m*(m+1)/2+j+1)%im)
          else
            call random_number(cc(m*(this%maxj+1)-m*(m+1)/2+j+1)%re)
          end if
        end do
      end do
    
    allocate( pmm(step), pmj(step), pmj1(step), pmj2(step), cosx(step), weight(step), sinx(step), ssym(step), asym(step), &
            & sumN(step*(1+this%maxj)), sumS(step*(1+this%maxj))                                                          )
      
      do
        cr = czero
        
        do i = 1, this%nLegendre, step
          call lege_setup_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
          call zero_poly_sub( step*(this%maxj+1), sumN(1), sumS(1) )
          
          do m = 0, this%maxj
            call pmm_recursion_sub( m, sinx(1), pmm(1) ) ; if (maxval(abs(pmm)) < this%tolm) exit
            
            call zero_poly_sub( step, ssym(1), asym(1) )
            
            j = m
              mj = m*(this%maxj+1)-m*(m+1)/2+j
              
              call pmj_setup_sub( pmm(1), pmj2(1), pmj1(1), pmj(1) )
              call pmj_backward_sub( 1, pmj(1), cc(1+mj), ssym(1) )
            
            do j = 1, (this%maxj-m)/2
              mj = mj+2
              
              call pmj_recursion_sub( this%amjrr(mj), this%bmjrr(mj), cosx(1), pmj2(1), pmj1(1), pmj(1) )
              call pmj_backward_sub( 1, pmj(1), cc(mj), asym(1) )
              
              call pmj_recursion_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
              call pmj_backward_sub( 1, pmj(1), cc(1+mj), ssym(1) )
              
              if ( maxval(abs(pmj)) < this%tolm ) exit
            end do
            
            if ( (maxval(abs(pmj)) >= this%tolm) .and. (mod((this%maxj-m),2) /= 0) ) then
              mj = mj+1
              
              call pmj_recursion_sub( this%amjrr(mj+1), this%bmjrr(mj+1), cosx(1), pmj2(1), pmj1(1), pmj(1) )
              call pmj_backward_sub( 1, pmj(1), cc(1+mj), asym(1) )
            end if
            
            call pmj_backward_recomb_sub( m, 1, ssym(1), asym(1), sumN(1+step*m), sumS(1+step*m) )
          end do
          
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

        maxdiff = 0._dbl
        
        do m = 0, this%jmax
          do j = m, this%jmax
            diff = abs( abs( cc(m*(this%maxj+1)-m*(m+1)/2+j+1) / cr(m*(this%maxj)-m*(m+1)/2+j+1) ) / this%scale / s4pi - 1 )
            
            if ( diff > maxdiff ) maxdiff = diff
          end do
        end do
        
        if ( maxdiff < 1.0d-5 ) then
          exit
        else
          this%tolm = this%tolm / 10
        end if
      end do
      
    deallocate( cr, cc, sumN, sumS, pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym )
    
  end subroutine vctol_sub
  
end submodule vctol