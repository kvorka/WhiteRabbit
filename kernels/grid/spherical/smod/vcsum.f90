submodule (SphericalHarmonics) vcsum
  implicit none ; contains
  
  module pure subroutine vcsum_sub(this, cajm, cbjm, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: cc(:), cr(:)
    
    !Array preparation
    allocate( cc(2*this%jms2) ); call zero_carray_sub( 2*this%jms2, cc(1) )
    allocate( cr(this%jms2)   ); call zero_carray_sub(   this%jms2, cr(1) )
    
    call this%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 2, 1 )
    call this%scal2scal_jm_to_mj_sub( cbjm(1), cc(1), 2, 2 )
    
    !Transform
    call this%lege_transform_sub( 1, 2, cc(1), cr(1), grid_op_2_vcsum_sub, grid_op_4_vcsum_sub, &
                                                    & grid_op_8_vcsum_sub, grid_op_16_vcsum_sub )
    
    !Rearranging indexing
    call this%scal2scal_mj_to_jm_sub( cr(1), 1, 1, cjm(1), 1, 1)
    
    !Cleaning
    deallocate( cr, cc )
    
    !Rescaling
    call this%rescale_sub( cjm(1), this%jms )
    
  end subroutine vcsum_sub
  
end submodule vcsum