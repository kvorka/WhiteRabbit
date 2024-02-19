submodule (SphericalHarmonics) vcsum
  implicit none ; contains
  
  module pure subroutine vcsum_sub(this, cajm, cbjm, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    integer                           :: i, j, m
    real(kind=dbl),       allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:), sinx(:), grid(:)
    complex(kind=dbl),    allocatable :: cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    allocate( cc(2*this%jms2), cr(this%jms2) )
      
      cc = czero
      cr = czero
      
      do m = 0, this%jmax
        do j = m, this%jmax
          cc(1+2*(m*(this%maxj+1)-m*(m+1)/2+j)) = cajm(j*(j+1)/2+m+1)
          cc(2+2*(m*(this%maxj+1)-m*(m+1)/2+j)) = cbjm(j*(j+1)/2+m+1)
        end do
      end do
      
    !Stepping of the algorithm :: 16
    allocate( pmm(16), pmj(16), pmj1(16), pmj2(16), cosx(16), weight(16), sinx(16), ssym(32), asym(32), &
            & sumN(32*(1+this%maxj)), sumS(32*(1+this%maxj)), grid(32*this%nFourier)                    )
      
      do i = 1, (this%nLegendre/16)*16, 16
        call lege_setup_16_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 32*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_16_sub( 2, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                         & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_16_vcsum_sub( grid(1), sumN(1) )
        call this%grid_op_16_vcsum_sub( grid(1), sumS(1) )
        
        call this%partial_forward_16_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                        & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
      !Stepping of the algorithm :: 8
      do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
        call lege_setup_8_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 16*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_8_sub( 2, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_8_vcsum_sub( grid(1), sumN(1) )
        call this%grid_op_8_vcsum_sub( grid(1), sumS(1) )
        
        call this%partial_forward_8_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
      !Stepping of the algorithm :: 4
      do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
        call lege_setup_4_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 8*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_4_sub( 2, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_4_vcsum_sub( grid(1), sumN(1) )
        call this%grid_op_4_vcsum_sub( grid(1), sumS(1) )
        
        call this%partial_forward_4_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
      !Stepping of the algorithm :: 2
      do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
        call lege_setup_2_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 4*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_2_sub( 2, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_2_vcsum_sub( grid(1), sumN(1) )
        call this%grid_op_2_vcsum_sub( grid(1), sumS(1) )
        
        call this%partial_forward_2_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
    deallocate( cc, sumN, sumS, grid, pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym )
      
      do j = 0, this%jmax
        m = 0
          cjm(j*(j+1)/2+m+1)%re = cr(m*(this%maxj+1)-m*(m+1)/2+j+1)%re * this%scale
          cjm(j*(j+1)/2+m+1)%im = 0._dbl
        
        do m = 1, j
          cjm(j*(j+1)/2+m+1) = cr(m*(this%maxj+1)-m*(m+1)/2+j+1) * this%scale
        end do
      end do
      
    deallocate(cr)
    
  end subroutine vcsum_sub
  
end submodule vcsum