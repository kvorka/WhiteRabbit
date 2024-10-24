submodule (lateral_grid) vcss_add_vcvv
  implicit none; contains
  
  module pure subroutine vcss_add_vcvv_sub(this, cajm, cbjm, cajml, cbjml, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*), cajml(*), cbjml(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_vectors_sub( 2, ca )
    call this%reindexing%allocate_vectors_sub( 8, cc )
    call this%reindexing%allocate_vectors_sub( 1, cr )
    
    call this%reindexing%vec2vec_jml_to_jml_sub( cajml(1), ca(1), 2, 1 )
    call this%reindexing%vec2vec_jml_to_jml_sub( cbjml(1), ca(1), 2, 2 )
    
    call this%reindexing%vec2scal_jml_to_mj_sub( ca(1), 2, cc(1), 8, 1 )
    call this%reindexing%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 8, 7 )
    call this%reindexing%scal2scal_jm_to_mj_sub( cbjm(1), cc(1), 8, 8 )
    
    deallocate(ca)
    
    !Transform
    call this%transform_sub( 1, 8, cc(1), cr(1), grid_op_vcss_add_vcvv_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2scal_mj_to_jm_sub( cr(1), 1, 1, cjm(1), 1, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end subroutine vcss_add_vcvv_sub
  
end submodule vcss_add_vcvv