submodule (SphericalHarmonics) vcsum
  implicit none ; contains
  
  pure subroutine grid_op_vcsum_sub(nfour, grid)
    integer,                intent(in)    :: nfour
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    gin(1:2,1:4,1:nfour) => grid(1:8*nfour)
    gout(1:4,1:nfour)    => grid(1:4*nfour)
    
    do i = 1, nfour
      do i2 = 1, 4
        gout(i2,i) = gin(1,i2,i) * gin(2,i2,i)
      end do
    end do
    
  end subroutine grid_op_vcsum_sub
  
  module subroutine vcsum_sub(this, cajm, cbjm, cjm)
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
    call this%lege_transform_sub( 1, 2, cc(1), cr(1), grid_op_vcsum_sub )
    
    !Rearranging indexing
    call this%scal2scal_mj_to_jm_sub( cr(1), 1, 1, cjm(1), 1, 1)
    
    !Cleaning
    deallocate( cr, cc )
    
  end subroutine vcsum_sub
  
end submodule vcsum