submodule (SphericalHarmonics) vcvv
  implicit none ; contains
  
  module pure subroutine vcvv_sub(this, cajml, cbjml, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajml(*), cbjml(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    integer                           :: i
    real(kind=dbl),       allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), weight(:), sinx(:), grid(:)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    !Array preparation
    allocate( ca(2*this%jmv1) ); call zero_carray_sub( 2*this%jmv1, ca(1) )
    allocate( cc(6*this%jms2) ); call zero_carray_sub( 6*this%jms2, cc(1) )
    allocate( cr(  this%jms2) ); call zero_carray_sub(   this%jms2, cr(1) )
      
      call this%vec2vec_jml_to_jml_sub( cajml(1), ca(1), 2, 1 )
      call this%vec2vec_jml_to_jml_sub( cbjml(1), ca(1), 2, 2 )
      
      call this%vec2scal_jml_to_mj_sub( ca(1), 2, cc(1) )
      
    deallocate(ca)
    
    !Allocating needed memory :: no reallocate for lower stepping
    allocate( pmm(16), pmj(16), pmj1(16), pmj2(16), cosx(16), weight(16), sinx(16), ssym(96),  &
            & asym(96), sumN(96*this%jmax3), sumS(96*this%jmax3), grid(96*this%nFourier)       )
      
      !Stepping of the algorithm :: 16
      do i = 1, (this%nLegendre/16)*16, 16
        call zero_carray_sub( 96*this%jmax3, sumN(1) )
        call zero_carray_sub( 96*this%jmax3, sumS(1) )

        call this%lege_init_16_sub( i-1, cosx(1), sinx(1), weight(1) )
        
        call this%partial_backward_16_sub( 6, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                         & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_16_vcvv_sub( grid(1), sumN(1) )
        call this%grid_op_16_vcvv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_16_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                        & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
      !Stepping of the algorithm :: 8
      do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
        call zero_carray_sub( 48*this%jmax3, sumN(1) )
        call zero_carray_sub( 48*this%jmax3, sumS(1) )
        
        call this%lege_init_8_sub( i-1, cosx(1), sinx(1), weight(1) )
        
        call this%partial_backward_8_sub( 6, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_8_vcvv_sub( grid(1), sumN(1) )
        call this%grid_op_8_vcvv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_8_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
      !Stepping of the algorithm :: 4
      do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
        call zero_carray_sub( 24*this%jmax3, sumN(1) )
        call zero_carray_sub( 24*this%jmax3, sumS(1) )

        call this%lege_init_4_sub( i-1, cosx(1), sinx(1), weight(1) )
        
        call this%partial_backward_4_sub( 6, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_4_vcvv_sub( grid(1), sumN(1) )
        call this%grid_op_4_vcvv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_4_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
      !Stepping of the algorithm :: 2
      do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
        call zero_carray_sub( 12*this%jmax3, sumN(1) )
        call zero_carray_sub( 12*this%jmax3, sumS(1) )
        
        call this%lege_init_2_sub( i-1, cosx(1), sinx(1), weight(1) )
        
        call this%partial_backward_2_sub( 6, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_2_vcvv_sub( grid(1), sumN(1) )
        call this%grid_op_2_vcvv_sub( grid(1), sumS(1) )
        
        call this%partial_forward_2_sub( 1, weight(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
    deallocate( cc, sumN, sumS, grid, pmm, pmj, pmj1, pmj2, cosx, sinx, weight, ssym, asym )
      
      !Rearranging indexing
      call this%scal2scal_mj_to_jm_sub( cr(1), 1, 1, cjm(1), 1, 1 )
      
    deallocate(cr)
    
    !Rescaling
    call this%rescale_sub( cjm(1), this%jms )
    
  end subroutine vcvv_sub
  
end submodule vcvv