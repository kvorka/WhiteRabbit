submodule (lateral_grid) vcvv
  implicit none ; contains
  
  module procedure vcvv_sub
    complex(kind=dbl), allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_vectors_sub( 2, ca )
    call this%reindexing%allocate_scalars_sub( 6, cc )
    call this%reindexing%allocate_scalars_sub( 1, cr )
    
    call this%reindexing%vec2vec_jml_to_jml_sub( cajml, ca, 2, 1 )
    call this%reindexing%vec2vec_jml_to_jml_sub( cbjml, ca, 2, 2 )
    call this%reindexing%vec2scal_jml_to_mj_sub( ca, 2, cc, 6, 1 )
    
    deallocate(ca)
    
    !Transform
    call this%transform_sub( 1, 6, cc, cr, grid_op_vcvv_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2scal_mj_to_jm_sub( cr, 1, 1, cjm, 1, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end procedure vcvv_sub
  
end submodule vcvv