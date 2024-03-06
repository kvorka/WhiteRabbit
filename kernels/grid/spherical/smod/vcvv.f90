submodule (SphericalHarmonics) vcvv
  implicit none ; contains
  
  module pure subroutine vcvv_sub(this, cajml, cbjml, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajml(*), cbjml(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    allocate( ca(2*this%jmv1) ); call zero_carray_sub( 2*this%jmv1, ca(1) )
    allocate( cc(6*this%jms2) ); call zero_carray_sub( 6*this%jms2, cc(1) )
    allocate( cr(  this%jms2) ); call zero_carray_sub(   this%jms2, cr(1) )
    
    call this%vec2vec_jml_to_jml_sub( cajml(1), ca(1), 2, 1 )
    call this%vec2vec_jml_to_jml_sub( cbjml(1), ca(1), 2, 2 )
    
    call this%vec2scal_jml_to_mj_sub( ca(1), 2, cc(1) )
    
    deallocate(ca)
    
    !Transform
    call this%lege_transform_sub( 1, 6, cc(1), cr(1), grid_op_2_vcvv_sub, grid_op_4_vcvv_sub, &
                                                    & grid_op_8_vcvv_sub, grid_op_16_vcvv_sub )
    
    !Rearranging indexing
    call this%scal2scal_mj_to_jm_sub( cr(1), 1, 1, cjm(1), 1, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
    !Rescaling
    call this%rescale_sub( cjm(1), this%jms )
    
  end subroutine vcvv_sub
  
end submodule vcvv