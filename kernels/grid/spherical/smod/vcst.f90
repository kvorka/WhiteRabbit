submodule (SphericalHarmonics) vcst
  implicit none ; contains
  
  module pure subroutine vcst_sub(this, cajm, cbjml2, cjml2)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjml2(*)
    complex(kind=dbl),    intent(out) :: cjml2(*)
    complex(kind=dbl),    allocatable :: cc(:), cr(:)
    
    !Array preparation
    allocate( cc(6*this%jms2) ); call zero_carray_sub( 6*this%jms2, cc(1) )
    allocate( cr(5*this%jms2) ); call zero_carray_sub( 5*this%jms2, cr(1) )
    
    call this%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 6, 1 )
    call this%devtens2scal_jml2_to_mj_sub( cbjml2(1), cc(1), 6, 2 )
    
    !Transform
    call this%lege_transform_sub( 5, 6, cc(1), cr(1), grid_op_2_vcst_sub, grid_op_4_vcst_sub, &
                                                    & grid_op_8_vcst_sub, grid_op_16_vcst_sub )
    
    !Rearranging indexing
    call this%scal2devtens_mj_to_jml2_sub( cr(1), 5, 1, cjml2(1) )
    
    !Cleaning
    deallocate( cc, cr )
    
    !Rescaling
    call this%rescale_sub( cjml2(1), this%jmt )
    
  end subroutine vcst_sub
  
end submodule vcst