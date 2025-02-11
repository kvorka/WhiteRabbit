submodule (lateral_grid) vcvxv
  implicit none; contains
  
  module procedure vcvxv_sub
    complex(kind=dbl), allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_vectors_sub( 2, ca )
    call this%reindexing%allocate_scalars_sub( 6, cc )
    call this%reindexing%allocate_scalars_sub( 3, cr )
    
    call this%reindexing%vec2vec_jml_to_jml_sub( cajml, ca, 2, 1 )
    call this%reindexing%vec2vec_jml_to_jml_sub( cbjml, ca, 2, 2 )
    call this%reindexing%vec2scal_jml_to_mj_sub( ca, 2, cc, 6, 1 )
    
    deallocate(ca)
    
    !Transform
    call this%transform_sub( 3, 6, cc, cr, grid_op_vcvxv_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2vec_mj_to_jml_sub( cr, 3, 1, cjml, 1, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end procedure vcvxv_sub
  
end submodule vcvxv