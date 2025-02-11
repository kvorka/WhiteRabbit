submodule (lateral_grid) vcsv
  implicit none; contains
  
  module procedure vcsv_sub
    complex(kind=dbl), allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_vectors_sub( 1, ca )
    call this%reindexing%allocate_scalars_sub( 4, cc )
    call this%reindexing%allocate_scalars_sub( 3, cr )
    
    call this%reindexing%scal2scal_jm_to_mj_sub( cajm, cc, 4, 1 )
    call this%reindexing%vec2vec_jml_to_jml_sub( cbjml, ca, 1, 1 )
    call this%reindexing%vec2scal_jml_to_mj_sub( ca, 1, cc, 4, 2 )
    
    deallocate(ca)
    
    !Transform
    call this%transform_sub( 3, 4, cc, cr, grid_op_vcsv_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2vecscal_mj_to_jm_sub( cr, 3, 1, cjm, 3, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end procedure vcsv_sub
  
end submodule vcsv