submodule (lateral_grid) vcvv_vcvgv
  implicit none; contains
  
  pure subroutine grid_op_vcvv_vcvgv_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2, i3
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:,:)
    real(kind=dbl), allocatable           :: tmp11(:), tmp12(:), tmp13(:), tmp(:)
    
    allocate( tmp11(step), tmp12(step), tmp13(step), tmp(step) )
    
    gin(1:step,1:3,1:5,1:nfour) => grid(1:15*step*nfour)
    gout(1:step,1:4,1:nfour)    => grid(1: 4*step*nfour)
    
    do i1 = 1, nfour
      tmp11 = gin(1:step,1,1,i1)
      tmp12 = gin(1:step,2,1,i1)
      tmp13 = gin(1:step,3,1,i1)
      
      do i3 = 1, 4
        tmp = gin(1:step,1,i3+1,i1)
          do concurrent ( i2 = 1:step )
            gout(i2,i3,i1) = tmp11(i2) * tmp(i2)
          end do
        
        tmp = gin(1:step,2,i3+1,i1)
          do concurrent ( i2 = 1:step )
            gout(i2,i3,i1) = gout(i2,i3,i1) + tmp12(i2) * tmp(i2)
          end do
        
        tmp = gin(1:step,3,i3+1,i1)
          do concurrent ( i2 = 1:step )
            gout(i2,i3,i1) = gout(i2,i3,i1) + tmp13(i2) * tmp(i2)
          end do
      end do
    end do
    
    deallocate( tmp11, tmp12, tmp13, tmp )
    
  end subroutine grid_op_vcvv_vcvgv_sub
  
  module pure subroutine vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(*), q(*), v(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_vectors_sub(  5, ca )
    call this%reindexing%allocate_scalars_sub( 15, cc )
    call this%reindexing%allocate_scalars_sub(  4, cr )
    
    call this%reindexing%vec2vec_jml_to_jml_sub( v(1), ca(1), 5, 1 )
    call this%reindexing%vec2vec_jml_to_jml_sub( q(1), ca(1), 5, 2 )
    call this%reindexing%gradvec2vec_jmlk_to_jml_sub( ri, v(1), dv_r(1), ca(1), 5, 3 )
    
    call this%reindexing%vec2scal_jml_to_mj_sub( ca(1), 5, cc(1), 15, 1 )
    
    deallocate(ca)
    
    !Transform
    call this%transform_sub( 4, 15, cc(1), cr(1), grid_op_vcvv_vcvgv_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2scal_mj_to_jm_sub( cr(1), 4, 1, cjm(1), 4, 1 )
    call this%reindexing%scal2vecscal_mj_to_jm_sub( cr(1), 4, 2, cjm(1), 4, 2 )
    
    !Cleaning
    deallocate( cc, cr )
    
  end subroutine vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv