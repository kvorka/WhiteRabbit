submodule (lateral_grid) vcvv_vcvgv
  implicit none; contains
  
  pure subroutine grid_op_vcvv_vcvgv_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:,:)
    real(kind=dbl), allocatable           :: tmp(:)
    
    allocate( tmp(3) )
    
    gin(1:3,1:5,1:step,1:nfour) => grid(1:15*step*nfour)
    gout(1:4,1:step,1:nfour)    => grid(1: 4*step*nfour)
    
    do i1 = 1, nfour
      do i2 = 1, step
        tmp = gin(1:3,1,i2,i1)
        
        gout(1,i2,i1) = sum( gin(1:3,2,i2,i1) * tmp(1:3) )
        gout(2,i2,i1) = sum( gin(1:3,3,i2,i1) * tmp(1:3) )
        gout(3,i2,i1) = sum( gin(1:3,4,i2,i1) * tmp(1:3) )
        gout(4,i2,i1) = sum( gin(1:3,5,i2,i1) * tmp(1:3) )
      end do
    end do
    
    deallocate( tmp )
    
  end subroutine grid_op_vcvv_vcvgv_sub
  
  module pure subroutine vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(*), q(*), v(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    allocate( ca( 5*this%reindexing%jmv1) ); call zero_carray_sub(  5*this%reindexing%jmv1, ca(1) ) 
    allocate( cc(15*this%reindexing%jms2) ); call zero_carray_sub( 15*this%reindexing%jms2, cc(1) ) 
    allocate( cr( 4*this%reindexing%jms2) ); call zero_carray_sub(  4*this%reindexing%jms2, cr(1) )
    
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