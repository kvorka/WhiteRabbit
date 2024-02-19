submodule (SphericalHarmonics) vctol
  implicit none ; contains
  
  module subroutine vctol_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer                             :: i, j, m, mj
    real(kind=dbl)                      :: diff, maxdiff, xrand
    real(kind=dbl),       allocatable   :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:), sinx(:)
    complex(kind=dbl),    allocatable   :: cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    allocate( cc(this%jms2), cr(this%jms2) )
      
      do m = 0, this%jmax+2
        do j = m, this%jmax+2
          mj = m*(this%maxj+1)-m*(m+1)/2+j+1
          
          if ( m /= 0 ) then
            call random_number( xrand ) ; cc(mj)%re = xrand
            call random_number( xrand ) ; cc(mj)%im = xrand
          else
            call random_number( xrand ) ; cc(mj)%re = xrand
            cc(mj)%im = czero
          end if
        end do
      end do
    
    allocate( pmm(16), pmj(16), pmj1(16), pmj2(16), cosx(16), weight(16), sinx(16), ssym(16), asym(16), &
            & sumN(16*(1+this%maxj)), sumS(16*(1+this%maxj))                                            )
      
      do
        cr = czero
        
        !Stepping of the algorithm :: 16
        do i = 1, (this%nLegendre/16)*16, 16
          call lege_setup_16_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
          call zero_poly_sub( 16*(this%maxj+1), sumN(1), sumS(1) )
          
          call this%partial_backward_16_sub( 1, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj(1), &
                                           & ssym(1), asym(1), cc(1), sumN(1), sumS(1)               )
          
          call this%partial_forward_16_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj(1), ssym(1), asym(1), &
                                          & cr(1), sumN(1), sumS(1) )
        end do
        
        !Stepping of the algorithm :: 8
        do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
          call lege_setup_8_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
          call zero_poly_sub( 8*(this%maxj+1), sumN(1), sumS(1) )
          
          call this%partial_backward_8_sub( 1, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj(1), &
                                           & ssym(1), asym(1), cc(1), sumN(1), sumS(1)               )
          
          call this%partial_forward_8_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj(1), ssym(1), asym(1), &
                                         & cr(1), sumN(1), sumS(1) )
        end do
        
        !Stepping of the algorithm :: 4
        do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
          call lege_setup_4_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
          call zero_poly_sub( 4*(this%maxj+1), sumN(1), sumS(1) )
          
          call this%partial_backward_4_sub( 1, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj(1), &
                                          & ssym(1), asym(1), cc(1), sumN(1), sumS(1)               )
          
          call this%partial_forward_4_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj(1), ssym(1), asym(1), &
                                         & cr(1), sumN(1), sumS(1) )
        end do
        
        !Stepping of the algorithm :: 2
        do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
          call lege_setup_2_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
          call zero_poly_sub( 2*(this%maxj+1), sumN(1), sumS(1) )
          
          call this%partial_backward_2_sub( 1, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj(1), &
                                          & ssym(1), asym(1), cc(1), sumN(1), sumS(1)               )
          
          call this%partial_forward_2_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj(1), ssym(1), asym(1), &
                                         & cr(1), sumN(1), sumS(1) )
        end do
        
        maxdiff = 0._dbl
        
        do mj = 1, this%jms2
          diff = abs( abs( cc(mj) / cr(mj) ) / this%scale / sqrt(4*pi) - 1 )
            
          if ( diff > maxdiff ) maxdiff = diff
        end do
        
        if ( maxdiff < 1.0d-5 ) then
          exit
        else
          this%tolm = this%tolm / 10
        end if
        
        if ( this%tolm < 1.0d-100) then
          write(*,*) 'Problem with precision setting, try another run (that usually fixes it).'
          stop
        end if
      end do
      
    deallocate( cr, cc, sumN, sumS, pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym )
    
  end subroutine vctol_sub
  
end submodule vctol