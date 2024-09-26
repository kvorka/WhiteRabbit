submodule (SphericalHarmonics) vcsv
  implicit none; contains
  
  pure subroutine grid_op_vcsv_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2
    real(kind=dbl)                        :: tmp
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:)
    
    gin(1:4,1:step,1:nfour)  => grid(1:4*step*nfour)
    gout(1:3,1:step,1:nfour) => grid(1:3*step*nfour)
    
    do i1 = 1, nfour
      do i2 = 1, step
        tmp = gin(1,i2,i1)
        
        gout(1,i2,i1) = gin(2,i2,i1) * tmp
        gout(2,i2,i1) = gin(3,i2,i1) * tmp
        gout(3,i2,i1) = gin(4,i2,i1) * tmp
      end do
    end do
    
  end subroutine grid_op_vcsv_sub
  
  module pure subroutine vcsv_sub(this, cajm, cbjml, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjml(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    allocate( ca(  this%reindexing%jmv1) ); call zero_carray_sub(   this%reindexing%jmv1, ca(1) ) 
    allocate( cc(4*this%reindexing%jms2) ); call zero_carray_sub( 4*this%reindexing%jms2, cc(1) ) 
    allocate( cr(3*this%reindexing%jms2) ); call zero_carray_sub( 3*this%reindexing%jms2, cr(1) )
    
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