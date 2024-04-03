submodule (SphericalHarmonics) vcst
  implicit none ; contains
  
  pure subroutine grid_op_vcst_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2
    real(kind=dbl)                        :: tmp
    real(kind=dbl), pointer               :: gin(:,:,:), gout(:,:,:)
    
    gin(1:6,1:step,1:nfour)  => grid(1:6*step*nfour)
    gout(1:5,1:step,1:nfour) => grid(1:5*step*nfour)
    
    do i1 = 1, nfour
      do i2 = 1, step
        tmp = gin(1,i2,i1)
        
        gout(1,i2,i1) = gin(2,i2,i1) * tmp
        gout(2,i2,i1) = gin(3,i2,i1) * tmp
        gout(3,i2,i1) = gin(4,i2,i1) * tmp
        gout(4,i2,i1) = gin(5,i2,i1) * tmp
        gout(5,i2,i1) = gin(6,i2,i1) * tmp
      end do
    end do
    
  end subroutine grid_op_vcst_sub
  
  module pure subroutine vcst_sub(this, cajm, cbjml2, cjml2)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjml2(*)
    complex(kind=dbl),    intent(out) :: cjml2(*)
    complex(kind=dbl),    allocatable :: cc(:), cr(:)
    
    !Array preparation
    allocate( cc(6*this%reindexing%jms2) ); call zero_carray_sub( 6*this%reindexing%jms2, cc(1) )
    allocate( cr(5*this%reindexing%jms2) ); call zero_carray_sub( 5*this%reindexing%jms2, cr(1) )
    
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