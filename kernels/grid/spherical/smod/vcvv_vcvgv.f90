submodule (SphericalHarmonics) vcvv_vcvgv
  implicit none ; contains
  
  module pure subroutine vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(*), q(*), v(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    allocate( ca( 5*this%jmv1) ); call zero_carray_sub(  5*this%jmv1, ca(1) ) 
    allocate( cc(15*this%jms2) ); call zero_carray_sub( 15*this%jms2, cc(1) ) 
    allocate( cr( 4*this%jms2) ); call zero_carray_sub(  4*this%jms2, cr(1) )
    
    call this%vec2vec_jml_to_jml_sub( q(1), ca(1), 5, 1 )
    call this%vec2vec_jml_to_jml_sub( v(1), ca(1), 5, 2 )
    call this%gradvec2vec_jmlk_to_jml_sub( ri, v(1), dv_r(1), ca(1), 5, 3 )
    
    call this%vec2scal_jml_to_mj_sub( ca(1), 5, cc(1) )
    
    deallocate(ca)
    
    !Transform
    call this%lege_transform_sub( 4, 15, cc(1), cr(1), grid_op_2_vcvv_vcvgv_sub, grid_op_4_vcvv_vcvgv_sub, &
                                                     & grid_op_8_vcvv_vcvgv_sub, grid_op_16_vcvv_vcvgv_sub )
    
    !Rearranging indexing
    call this%scal2scal_mj_to_jm_sub( cr(1), 4, 1, cjm(1), 4, 1 )
    call this%scal2vecscal_mj_to_jm_sub( cr(1), 4, 2, cjm(1), 4, 2 )
    
    !Cleaning
    deallocate( cc, cr )
    
    !Rescaling
    call this%rescale_sub( cjm(1), 4*this%jms )
    
  end subroutine vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv