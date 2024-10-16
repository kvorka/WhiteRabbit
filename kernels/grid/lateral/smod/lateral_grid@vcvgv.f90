submodule (lateral_grid) vcvgv
  implicit none; contains
  
  pure subroutine grid_op_vcvgv_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:,:)
    real(kind=dbl), allocatable           :: tmp(:)
    
    gin(1:3,1:4,1:step,1:nfour) => grid(1:12*step*nfour)
    gout(1:3,1:step,1:nfour)    => grid(1: 3*step*nfour)
    
    allocate( tmp(3) )
    
    do i1 = 1, nfour
      do i2 = 1, step
        tmp = gin(1:3,1,i2,i1)
        
        gout(1,i2,i1) = sum( gin(1:3,2,i2,i1) * tmp(1:3) )
        gout(2,i2,i1) = sum( gin(1:3,3,i2,i1) * tmp(1:3) )
        gout(3,i2,i1) = sum( gin(1:3,4,i2,i1) * tmp(1:3) )
      end do
    end do
    
    deallocate( tmp )
    
  end subroutine grid_op_vcvgv_sub
  
  module pure subroutine vcvgv_sub(this, ri, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(*), v(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_vectors_sub(  4, ca )
    call this%reindexing%allocate_scalars_sub( 12, cc )
    call this%reindexing%allocate_scalars_sub(  3, cr )
    
    call this%reindexing%vec2vec_jml_to_jml_sub( v(1), ca(1), 4, 1 )
    call this%reindexing%gradvec2vec_jmlk_to_jml_sub( ri, v(1), dv_r(1), ca(1), 4, 2 )
    call this%reindexing%vec2scal_jml_to_mj_sub( ca(1), 4, cc(1), 12, 1 )
    
    deallocate(ca)
    
    !Transform
    call this%transform_sub( 3, 12, cc(1), cr(1), grid_op_vcvgv_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2vecscal_mj_to_jm_sub( cr(1), 3, 1, cjm(1), 3, 1 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end subroutine vcvgv_sub
  
end submodule vcvgv