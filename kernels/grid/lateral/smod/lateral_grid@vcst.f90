submodule (lateral_grid) vcst
  implicit none ; contains
  
  module procedure vcst_sub
    complex(kind=dbl), allocatable :: cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_scalars_sub( 6, cc )
    call this%reindexing%allocate_scalars_sub( 5, cr )
    
    call this%reindexing%scal2scal_jm_to_mj_sub( cajm, cc, 6, 1 )
    call this%reindexing%devtens2scal_jml2_to_mj_sub( cbjml2, cc, 6, 2 )
    
    !Transform
    call this%transform_sub( 5, 6, cc, cr, grid_op_vcst_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2devtens_mj_to_jml2_sub( cr, 5, 1, cjml2 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end procedure vcst_sub
  
end submodule vcst