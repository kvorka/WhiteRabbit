submodule (lateral_grid) vcst
  implicit none ; contains
  
  module pure subroutine vcst_sub(this, cajm, cbjml2, cjml2)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjml2(*)
    complex(kind=dbl),    intent(out) :: cjml2(*)
    complex(kind=dbl),    allocatable :: cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_scalars_sub( 6, cc )
    call this%reindexing%allocate_scalars_sub( 5, cr )
    
    call this%reindexing%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 6, 1 )
    call this%reindexing%devtens2scal_jml2_to_mj_sub( cbjml2(1), cc(1), 6, 2 )
    
    !Transform
    call this%transform_sub( 5, 6, cc(1), cr(1), grid_op_vcst_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2devtens_mj_to_jml2_sub( cr(1), 5, 1, cjml2(1) )
    
    !Cleaning
    deallocate( cc, cr )
    
  end subroutine vcst_sub
  
end submodule vcst