submodule (lateral_grid) vcsv
  implicit none; contains
  
  pure subroutine grid_op_vcsv_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2, i3
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:)
    real(kind=dbl), allocatable           :: tmp(:), tmp1(:)
    
    gin(1:step,1:4,1:nfour)  => grid(1:4*step*nfour)
    gout(1:step,1:3,1:nfour) => grid(1:3*step*nfour)
    
    allocate( tmp(step), tmp1(step) )
    
    do i1 = 1, nfour
      tmp = gin(1:step,1,i1)
      
      do i3 = 1, 3
        tmp1 = gin(1:step,i3+1,i1)
        
        do concurrent ( i2 = 1:step )
          gout(i2,i3,i1) = tmp(i2) * tmp1(i2)
        end do
      end do
    end do
    
    deallocate( tmp, tmp1 )
    
  end subroutine grid_op_vcsv_sub
  
  module pure subroutine vcsv_sub(this, cajm, cbjml, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjml(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_vectors_sub( 1, ca )
    call this%reindexing%allocate_scalars_sub( 4, cc )
    call this%reindexing%allocate_scalars_sub( 3, cr )
    
    call this%reindexing%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 4, 1 )
    call this%reindexing%vec2vec_jml_to_jml_sub( cbjml(1), ca(1), 1, 1 )
    call this%reindexing%vec2scal_jml_to_mj_sub( ca(1), 1, cc(1), 4, 2 )
    
    deallocate(ca)
    
    !Transform
    call this%transform_sub( 3, 4, cc(1), cr(1), grid_op_vcsv_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2vecscal_mj_to_jm_sub( cr(1), 3, 1, cjm(1), 3, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end subroutine vcsv_sub
  
end submodule vcsv