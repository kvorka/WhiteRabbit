submodule (lateral_grid) vcss_add_vcvv
  implicit none; contains
  
  module procedure vcss_add_vcvv_sub
    complex(kind=dbl), allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_vectors_sub( 2, ca )
    call this%reindexing%allocate_vectors_sub( 8, cc )
    call this%reindexing%allocate_vectors_sub( 1, cr )
    
    call this%reindexing%vec2vec_jml_to_jml_sub( cajml, ca, 2, 1 )
    call this%reindexing%vec2vec_jml_to_jml_sub( cbjml, ca, 2, 2 )
    
    call this%reindexing%vec2scal_jml_to_mj_sub( ca, 2, cc, 8, 1 )
    call this%reindexing%scal2scal_jm_to_mj_sub( cajm, cc, 8, 7 )
    call this%reindexing%scal2scal_jm_to_mj_sub( cbjm, cc, 8, 8 )
    
    deallocate(ca)
    
    !Transform
    call this%transform_sub( 1, 8, cc, cr, grid_op_vcss_add_vcvv_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2scal_mj_to_jm_sub( cr, 1, 1, cjm, 1, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end procedure vcss_add_vcvv_sub
  
end submodule vcss_add_vcvv