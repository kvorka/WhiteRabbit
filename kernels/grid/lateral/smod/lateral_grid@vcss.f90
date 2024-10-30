submodule (lateral_grid) vcss
  implicit none ; contains
  
  module procedure vcss_sub
    complex(kind=dbl), allocatable :: cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_scalars_sub( 2, cc )
    call this%reindexing%allocate_scalars_sub( 1, cr )
    
    call this%reindexing%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 2, 1 )
    call this%reindexing%scal2scal_jm_to_mj_sub( cbjm(1), cc(1), 2, 2 )
    
    !Transform
    call this%transform_sub( 1, 2, cc(1), cr(1), grid_op_vcss_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2scal_mj_to_jm_sub( cr(1), 1, 1, cjm(1), 1, 1)
    
    !Cleaning
    deallocate( cr, cc )
    
  end procedure vcss_sub
  
end submodule vcss