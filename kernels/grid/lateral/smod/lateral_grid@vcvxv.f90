submodule (lateral_grid) vcvxv
  implicit none; contains
  
  pure subroutine grid_op_vcvxv_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:,:)
    real(kind=dbl), allocatable           :: tmp11(:), tmp12(:), tmp21(:), tmp22(:)
    
    gin(1:step,1:3,1:2,1:nfour) => grid(1:6*step*nfour)
    gout(1:step,1:3,1:nfour)    => grid(1:3*step*nfour)
    
    allocate( tmp11(step), tmp12(step), tmp21(step), tmp22(step) )
    
    do i1 = 1, nfour
      tmp11 = gin(1:step,2,1,i1)
      tmp12 = gin(1:step,3,1,i1)
      tmp21 = gin(1:step,2,2,i1)
      tmp22 = gin(1:step,3,2,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,1,i1) = tmp11(i2) * tmp22(i2) - tmp12(i2) * tmp21(i2)
      end do
      
      tmp11 = gin(1:step,3,1,i1)
      tmp12 = gin(1:step,1,1,i1)
      tmp21 = gin(1:step,3,2,i1)
      tmp22 = gin(1:step,1,2,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,2,i1) = tmp11(i2) * tmp22(i2) - tmp12(i2) * tmp21(i2)
      end do
      
      tmp11 = gin(1:step,1,1,i1)
      tmp12 = gin(1:step,2,1,i1)
      tmp21 = gin(1:step,1,2,i1)
      tmp22 = gin(1:step,2,2,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,3,i1) = tmp11(i2) * tmp22(i2) - tmp12(i2) * tmp21(i2)
      end do
    end do
    
    deallocate( tmp11, tmp12, tmp21, tmp22 )
    
  end subroutine grid_op_vcvxv_sub
  
  module pure subroutine vcvxv_sub(this, cajml, cbjml, cjml)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajml(*), cbjml(*)
    complex(kind=dbl),    intent(out) :: cjml(*)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_vectors_sub( 2, ca )
    call this%reindexing%allocate_scalars_sub( 6, cc )
    call this%reindexing%allocate_scalars_sub( 3, cr )
    
    call this%reindexing%vec2vec_jml_to_jml_sub( cajml(1), ca(1), 2, 1 )
    call this%reindexing%vec2vec_jml_to_jml_sub( cbjml(1), ca(1), 2, 2 )
    call this%reindexing%vec2scal_jml_to_mj_sub( ca(1), 2, cc(1), 6, 1 )
    
    deallocate(ca)
    
    !Transform
    call this%transform_sub( 3, 6, cc(1), cr(1), grid_op_vcvxv_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2vec_mj_to_jml_sub( cr(1), 3, 1, cjml(1), 1, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end subroutine vcvxv_sub
  
end submodule vcvxv