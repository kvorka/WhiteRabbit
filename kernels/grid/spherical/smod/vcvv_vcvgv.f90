submodule (SphericalHarmonics) vcvv_vcvgv
  implicit none ; contains
  
  module pure subroutine vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(*), q(*), v(*)
    complex(kind=dbl),    intent(out) :: cjm(4,*)
    integer                           :: i, j, m, ijm, mj, mj1, mj2
    real(kind=dbl),       allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), sinx(:), weight(:), grid(:)
    complex(kind=dbl),    allocatable :: cab(:,:), cc(:), cr(:,:), ssym(:), asym(:), sumN(:), sumS(:)
    
    !Preparing arrays for transforms
    allocate( cab(5,this%jmv1), cc(15*this%jms2), cr(4,this%jms2) )
      
      cab = czero
      cc  = czero
      cr  = czero
      
      do concurrent( ijm = 1:this%jmv )
        cab(1,ijm) = q(ijm)
        cab(2,ijm) = v(ijm)
      end do
      
      call gradvec_to_vectors_sub( this%jmax, 3, 5, ri, v(1), dv_r(1), cab(1,1) )
      call vectors_to_scalars_sub( this%jmax, 5, cab(1,1), cc(1) )
      
    deallocate(cab)

    !Allocating needed memory :: no reallocate for lower stepping
    allocate( pmm(16), pmj(16), pmj1(16), pmj2(16), cosx(16), weight(16), sinx(16), ssym(240),     &
            & asym(240), sumN(240*(this%maxj+1)), sumS(240*(this%maxj+1)), grid(240*this%nFourier) )
    
      !Stepping of the algorithm :: 16
      do i = 1, (this%nLegendre/16)*16, 16
        call lege_setup_16_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 240*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_16_sub( 15, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),   &
                                         & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_16_vcvv_vcvgv_sub( grid(1), sumN(1) )
        call this%grid_op_16_vcvv_vcvgv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_16_sub( 4, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                        & pmj(1), ssym(1), asym(1), cr(1,1), sumN(1), sumS(1)       )
      end do
    
      !Stepping of the algorithm :: 8
      do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
        call lege_setup_8_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 120*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_8_sub( 15, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),   &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_8_vcvv_vcvgv_sub( grid(1), sumN(1) )
        call this%grid_op_8_vcvv_vcvgv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_8_sub( 4, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1,1), sumN(1), sumS(1)       )
      end do
      
      !Stepping of the algorithm :: 4
      do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
        call lege_setup_4_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 60*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_4_sub( 15, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),   &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_4_vcvv_vcvgv_sub( grid(1), sumN(1) )
        call this%grid_op_4_vcvv_vcvgv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_4_sub( 4, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1,1), sumN(1), sumS(1)       )
      end do
      
      !Stepping of the algorithm :: 2
      do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
        call lege_setup_2_sub( this%roots(i), this%fftLege(i), cosx(1), sinx(1), weight(1) )
        call zero_poly_sub( 30*(this%maxj+1), sumN(1), sumS(1) )
        
        call this%partial_backward_2_sub( 15, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),   &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_2_vcvv_vcvgv_sub( grid(1), sumN(1) )
        call this%grid_op_2_vcvv_vcvgv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_2_sub( 4, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1,1), sumN(1), sumS(1)       )
      end do
    
    deallocate( cc, sumN, sumS, grid, pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym )
      
      call cartesian_to_cyclic_sub( 2, 4, this%jms2, cr(1,1) )
      
      j = 0
        m = 0
        ijm = 1
          mj  = m*(this%maxj+1)-m*(m+1)/2+j+1
          mj2 = mj + this%maxj + 1 - m
          
          cjm(1,ijm) =        cr(1,mj  )                                ; cjm(1,ijm)%im = zero
          cjm(4,ijm) =        cr(2,mj2 )   * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                     &        cr(4,mj+1)   * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                     & conjg( cr(2,mj2 ) ) * cleb1_fn(j+1,m-1,1,+1,j,m) ; cjm(4,ijm)%im = zero
    
      do j = 1, this%jmax
        m = 0
          ijm = ijm+1
          mj  = m*(this%maxj+1)-m*(m+1)/2+j+1
          mj2 = mj + this%maxj - m
          
          cjm(1,ijm) =        cr(1,mj   )                                ; cjm(1,ijm)%im = zero
          cjm(2,ijm) =        cr(2,mj2-1)   * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                     &        cr(4,mj -1)   * cleb1_fn(j-1,m+0,1, 0,j,m) + &
                     & conjg( cr(2,mj2-1) ) * cleb1_fn(j-1,m-1,1,+1,j,m) ; cjm(2,ijm)%im = zero
          cjm(3,ijm) =        cr(2,mj2  )   * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                     &        cr(4,mj   )   * cleb1_fn(j  ,m+0,1, 0,j,m) + &
                     & conjg( cr(2,mj2  ) ) * cleb1_fn(j  ,m-1,1,+1,j,m) ; cjm(3,ijm)%re = zero
          cjm(4,ijm) =        cr(2,mj2+1)   * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                     &        cr(4,mj +1)   * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                     & conjg( cr(2,mj2+1) ) * cleb1_fn(j+1,m-1,1,+1,j,m) ; cjm(4,ijm)%im = zero
        
        do m = 1, j
          ijm = ijm+1
          mj  = m*(this%maxj+1)-m*(m+1)/2+j+1
          mj1 = mj - this%maxj - 1 + m
          mj2 = mj + this%maxj - m
          
          cjm(1,ijm) = cr(1,mj   )
          cjm(2,ijm) = cr(2,mj2-1) * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                     & cr(4,mj -1) * cleb1_fn(j-1,m+0,1, 0,j,m) + &
                     & cr(3,mj1-1) * cleb1_fn(j-1,m-1,1,+1,j,m)
          cjm(3,ijm) = cr(2,mj2  ) * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                     & cr(4,mj   ) * cleb1_fn(j  ,m+0,1, 0,j,m) + &
                     & cr(3,mj1  ) * cleb1_fn(j  ,m-1,1,+1,j,m)
          cjm(4,ijm) = cr(2,mj2+1) * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                     & cr(4,mj +1) * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                     & cr(3,mj1+1) * cleb1_fn(j+1,m-1,1,+1,j,m)
        end do
      end do
      
    deallocate(cr)
    
    do concurrent ( ijm = 1:this%jms )
      cjm(1,ijm) = cjm(1,ijm) * this%scale
      cjm(2,ijm) = cjm(2,ijm) * this%scale
      cjm(3,ijm) = cjm(3,ijm) * this%scale
      cjm(4,ijm) = cjm(4,ijm) * this%scale
    end do
    
  end subroutine vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv