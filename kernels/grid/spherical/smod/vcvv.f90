submodule (SphericalHarmonics) vcvv
  implicit none ; contains
  
  pure subroutine grid_op_vcvv_sub(nfour, grid)
    integer,                intent(in)    :: nfour
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    gin(1:6,1:4,1:nfour) => grid(1:24*nfour)
    gout(1:4,1:nfour)    => grid(1: 4*nfour)
    
    do i = 1, nfour
      do i2 = 1, 4
        gout(i2,i) = gin(1,i2,i) * gin(4,i2,i) + gin(2,i2,i) * gin(5,i2,i) + gin(3,i2,i) * gin(6,i2,i)
      end do
    end do
    
  end subroutine grid_op_vcvv_sub
  
  module subroutine vcvv_sub(this, cajml, cbjml, cjm)
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
    call this%lege_transform_sub( 1, 6, cc(1), cr(1), grid_op_vcvv_sub )
    
    !Rearranging indexing
    call this%scal2scal_mj_to_jm_sub( cr(1), 1, 1, cjm(1), 1, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end subroutine vcvv_sub
  
end submodule vcvv