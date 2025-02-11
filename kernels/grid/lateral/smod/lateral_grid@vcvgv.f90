submodule (lateral_grid) vcvgv
  implicit none; contains
  
  module procedure vcvgv_sub
    complex(kind=dbl), allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_vectors_sub(  4, ca )
    call this%reindexing%allocate_scalars_sub( 12, cc )
    call this%reindexing%allocate_scalars_sub(  3, cr )
    
    call this%reindexing%vec2vec_jml_to_jml_sub( v, ca, 4, 1 )
    call this%reindexing%gradvec2vec_jmlk_to_jml_sub( ri, v, dv_r, ca, 4, 2 )
    call this%reindexing%vec2scal_jml_to_mj_sub( ca, 4, cc, 12, 1 )
    
    deallocate(ca)
    
    !Transform
    call this%transform_sub( 3, 12, cc, cr, grid_op_vcvgv_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2vecscal_mj_to_jm_sub( cr, 3, 1, cjm, 3, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end procedure vcvgv_sub
  
end submodule vcvgv