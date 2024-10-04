submodule (lateral_grid) vcss
  implicit none ; contains
  
  pure subroutine grid_op_vcss_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    gin(1:2,1:step,1:nfour) => grid(1:2*step*nfour)
    gout(1:step,1:nfour)    => grid(1:  step*nfour)
    
    do i1 = 1, nfour
      do i2 = 1, step
        gout(i2,i1) = gin(1,i2,i1) * gin(2,i2,i1)
      end do
    end do
    
  end subroutine grid_op_vcss_sub
  
  module pure subroutine vcss_sub(this, cajm, cbjm, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: cc(:), cr(:)
    
    !Array preparation
    allocate( cc(2*this%reindexing%jms2) ); call zero_carray_sub( 2*this%reindexing%jms2, cc(1) )
    allocate( cr(  this%reindexing%jms2) ); call zero_carray_sub(   this%reindexing%jms2, cr(1) )
    
    call this%reindexing%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 2, 1 )
    call this%reindexing%scal2scal_jm_to_mj_sub( cbjm(1), cc(1), 2, 2 )
    
    !Transform
    call this%transform_sub( 1, 2, cc(1), cr(1), grid_op_vcss_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2scal_mj_to_jm_sub( cr(1), 1, 1, cjm(1), 1, 1)
    
    !Cleaning
    deallocate( cr, cc )
    
  end subroutine vcss_sub
  
end submodule vcss