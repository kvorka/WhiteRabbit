submodule (lateral_grid) vcvxv
  implicit none; contains
  
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