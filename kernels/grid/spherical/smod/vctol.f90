submodule (SphericalHarmonics) vctol
  implicit none ; contains
  
  module subroutine vctol_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer                             :: i, j, m, mj
    real(kind=dbl)                      :: diff, maxdiff, xrand
    real(kind=dbl),       allocatable   :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), w(:), sinx(:)
    complex(kind=dbl),    allocatable   :: cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    allocate( cc(this%jms2), cr(this%jms2) )
      
      do m = 0, this%jmax2
        do j = m, this%jmax2
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          if ( m /= 0 ) then
            call random_number( xrand ) ; cc(mj)%re = xrand
            call random_number( xrand ) ; cc(mj)%im = xrand
          else
            call random_number( xrand ) ; cc(mj)%re = xrand
            cc(mj)%im = czero
          end if
        end do
      end do
    
    allocate( pmm(16), pmj(16), pmj1(16), pmj2(16), cosx(16), w(16), sinx(16), &
            & ssym(16), asym(16), sumN(16*this%jmax3), sumS(16*this%jmax3)          )
      
    do
      cr = czero
      
      !Stepping of the algorithm :: 16
      do i = 1, (this%nLegendre/16)*16, 16
        call zero_carray_sub( 16*this%jmax3, sumN(1) )
        call zero_carray_sub( 16*this%jmax3, sumS(1) )
        
        call this%lege_init_16_sub( i-1, cosx(1), sinx(1), w(1) )
        
        call this%partial_backward_16_sub( 1, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), pmj(1), &
                                         & ssym(1), asym(1), cc(1), sumN(1), sumS(1)               )
        
        call this%partial_forward_16_sub( 1, w(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                        & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
      !Stepping of the algorithm :: 8
      do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
        call zero_carray_sub( 8*this%jmax3, sumN(1) )
        call zero_carray_sub( 8*this%jmax3, sumS(1) )
        
        call this%lege_init_8_sub( i-1, cosx(1), sinx(1), w(1) )
        
        call this%partial_backward_8_sub( 1, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%partial_forward_8_sub( 1, w(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
      !Stepping of the algorithm :: 4
      do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
        call zero_carray_sub( 4*this%jmax3, sumN(1) )
        call zero_carray_sub( 4*this%jmax3, sumS(1) )
        
        call this%lege_init_4_sub( i-1, cosx(1), sinx(1), w(1) )
        
        call this%partial_backward_4_sub( i, 1, cosx(1), pmj2(1), pmj1(1), pmj(1),  &
                                        & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%partial_forward_4_sub( i, 1, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                       & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
      !Stepping of the algorithm :: 2
      do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
        call zero_carray_sub( 2*this%jmax3, sumN(1) )
        call zero_carray_sub( 2*this%jmax3, sumS(1) )
        
        call this%lege_init_2_sub( i-1, cosx(1), sinx(1), w(1) )
        
        call this%partial_backward_2_sub( i, 1, cosx(1), pmj2(1), pmj1(1), pmj(1),  &
                                        & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%partial_forward_2_sub( i, 1, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                       & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
      maxdiff = 0._dbl
      
      do mj = 1, this%jms2
        diff = abs( abs( cc(mj) / cr(mj) ) / this%scale / sqrt(4*pi) - 1 )
          
        if ( diff > maxdiff ) maxdiff = diff
      end do
      
      if ( maxdiff <= 1.0d-4 ) then
        exit
      else
        this%tolm = this%tolm / 10
      end if
      
      if ( this%tolm < 1.0d-100) then
        write(*,*) 'Problem with precision setting, try another run (that usually fixes it).'
        stop
      end if
    end do
      
    deallocate( cr, cc, sumN, sumS, pmm, pmj, pmj1, pmj2, cosx, sinx, w, ssym, asym )
    
  end subroutine vctol_sub
  
end submodule vctol
