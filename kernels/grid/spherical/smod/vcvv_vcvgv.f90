submodule (SphericalHarmonics) vcvv_vcvgv
  implicit none ; contains
  
  module pure subroutine vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(*), q(*), v(*)
    complex(kind=dbl),    intent(out) :: cjm(4,*)
    integer                           :: i
    real(kind=dbl),       allocatable :: pmj(:), pmj1(:), pmj2(:), cosx(:), w(:), grid(:)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    !Array preparation
    allocate( ca( 5*this%jmv1) ); call zero_carray_sub(  5*this%jmv1, ca(1) ) 
    allocate( cc(15*this%jms2) ); call zero_carray_sub( 15*this%jms2, cc(1) ) 
    allocate( cr( 4*this%jms2) ); call zero_carray_sub(  4*this%jms2, cr(1) )
      
      call this%vec2vec_jml_to_jml_sub( q(1), ca(1), 5, 1 )
      call this%vec2vec_jml_to_jml_sub( v(1), ca(1), 5, 2 )
      call this%gradvec2vec_jmlk_to_jml_sub( ri, v(1), dv_r(1), ca(1), 5, 3 )
      
      call this%vec2scal_jml_to_mj_sub( ca(1), 5, cc(1) )
      
    deallocate(ca)
    
    !Allocating needed memory :: no reallocate for lower stepping
    allocate( pmj(16), pmj1(16), pmj2(16), cosx(16), w(16), ssym(240), asym(240), &
            & sumN(240*this%jmax3), sumS(240*this%jmax3), grid(240*this%nFourier) )
      
      !Stepping of the algorithm :: 16
      do i = 1, (this%nLegendre/16)*16, 16
        call this%partial_backward_16_sub( i, 15, cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                         & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_16_vcvv_vcvgv_sub( grid(1), sumN(1) )
        call this%grid_op_16_vcvv_vcvgv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_16_sub( i, 4, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                        & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
      !Stepping of the algorithm :: 8
      do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
        call this%partial_backward_8_sub( i, 15, cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                        & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_8_vcvv_vcvgv_sub( grid(1), sumN(1) )
        call this%grid_op_8_vcvv_vcvgv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_8_sub( i, 4, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                       & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
      !Stepping of the algorithm :: 4
      do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
        call this%partial_backward_4_sub( i, 15, cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                        & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_4_vcvv_vcvgv_sub( grid(1), sumN(1) )
        call this%grid_op_4_vcvv_vcvgv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_4_sub( i, 4, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                       & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
      !Stepping of the algorithm :: 2
      do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
        call this%partial_backward_2_sub( i, 15, cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                        & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_2_vcvv_vcvgv_sub( grid(1), sumN(1) )
        call this%grid_op_2_vcvv_vcvgv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_2_sub( i, 4, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                       & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
    
    deallocate( cc, sumN, sumS, grid, pmj, pmj1, pmj2, cosx, w, ssym, asym )
      
      !Rearranging indexing
      call this%scal2scal_mj_to_jm_sub( cr(1), 4, 1, cjm(1,1), 4, 1 )
      call this%scal2vecscal_mj_to_jm_sub( cr(1), 4, 2, cjm(1,1), 4, 2 )
      
    deallocate(cr)
    
    !Rescaling
    call this%rescale_sub( cjm(1,1), 4*this%jms )
    
  end subroutine vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv