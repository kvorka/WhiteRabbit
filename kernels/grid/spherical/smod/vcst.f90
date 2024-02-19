submodule (SphericalHarmonics) vcst
  implicit none ; contains
  
  module pure subroutine vcst_sub(this, cajm, cbjml2, cjml2)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjml2(*)
    complex(kind=dbl),    intent(out) :: cjml2(*)
    integer                           :: i, j, m, l, k, lm
    complex(kind=dbl)                 :: csum
    real(kind=dbl),       allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:), sinx(:), grid(:)
    complex(kind=dbl),    allocatable :: ssym(:), asym(:), sumN(:), sumS(:), caux(:), cc(:), cr(:,:)
    
    allocate( cc(6*this%jms2), cr(5,this%jms2) )
      
      cc  = czero
      cr  = czero
      
      do j = 0, this%jmax
        do m = 0, j
          cc(6+6*(m*(this%jmax+3)-m*(m+1)/2+j)) = cajm(j*(j+1)/2+m+1)
        end do
      end do
      
      call tensor2_to_scalars_sub(this%jmax, 1, 6, cbjml2(1), cc(1))
    
    allocate( pmm(16), pmj(16), pmj1(16), pmj2(16), cosx(16), weight(16), sinx(16), ssym(96), asym(96), &
            & sumN(96*(1+this%maxj)), sumS(96*(1+this%maxj)), grid(96*this%nFourier)                    )
      
      !Stepping of the algorithm :: 16
      do i = 1, (this%nLegendre/16)*16, 16
        call lege_setup_16_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 96*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_16_sub( 6, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                         & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_16_vcst_sub( grid(1), sumN(1) )
        call this%grid_op_16_vcst_sub( grid(1), sumS(1) )
        
        call this%partial_forward_16_sub( 5, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                        & pmj(1), ssym(1), asym(1), cr(1,1), sumN(1), sumS(1)       )
      end do
      
      !Stepping of the algorithm :: 8
      do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
        call lege_setup_8_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 48*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_8_sub( 6, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_8_vcst_sub( grid(1), sumN(1) )
        call this%grid_op_8_vcst_sub( grid(1), sumS(1) )
        
        call this%partial_forward_8_sub( 5, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1,1), sumN(1), sumS(1)       )
      end do
      
      !Stepping of the algorithm :: 4
      do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
        call lege_setup_4_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 24*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_4_sub( 6, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_4_vcst_sub( grid(1), sumN(1) )
        call this%grid_op_4_vcst_sub( grid(1), sumS(1) )
        
        call this%partial_forward_4_sub( 5, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1,1), sumN(1), sumS(1)       )
      end do
      
      !Stepping of the algorithm :: 2
      do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
        call lege_setup_2_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 12*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_2_sub( 6, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_2_vcst_sub( grid(1), sumN(1) )
        call this%grid_op_2_vcst_sub( grid(1), sumS(1) )
        
        call this%partial_forward_2_sub( 5, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1,1), sumN(1), sumS(1)       )
      end do
      
    deallocate( cc, sumN, sumS, grid, pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym )
    allocate( caux(5) )
      
      do j = 0, this%jmax
        do m = 0, j
          do l = abs(j-2), j+2
            csum = czero
            
            do k = -2, 2
              lm = abs(m-k)*(this%maxj+1)-abs(m-k)*(abs(m-k)+1)/2+l+1
              
              if (m < k) then
                caux = (-1)**(k-m) * conjg( cr(:,lm) )
              else
                caux = cr(:,lm)
              end if
              
              select case (k)
                case (-2)
                  csum = csum + ( caux(1) + cunit * caux(5) ) * cleb2_fn(l,m+2,2,-2,j,m)
                case (-1)
                  csum = csum + ( caux(2) - cunit * caux(4) ) * cleb2_fn(l,m+1,2,-1,j,m)
                case ( 0)
                  csum = csum + ( caux(3) +         caux(3) ) * cleb2_fn(l,m  ,2, 0,j,m)
                case (+1)
                  csum = csum - ( caux(2) + cunit * caux(4) ) * cleb2_fn(l,m-1,2,+1,j,m)
                case (+2)
                  csum = csum + ( caux(1) - cunit * caux(5) ) * cleb2_fn(l,m-2,2,+2,j,m)
              end select
            end do
            
            cjml2(5*(j*(j+1)/2+m)+l-j-1) = this%scale * csum / 2
          end do
        end do
      end do
      
    deallocate( cr, caux )
    
  end subroutine vcst_sub
  
end submodule vcst