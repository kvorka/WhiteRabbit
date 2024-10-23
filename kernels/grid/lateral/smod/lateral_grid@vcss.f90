submodule (lateral_grid) vcss
  implicit none ; contains
  
  pure subroutine grid_op_vcss_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    real(kind=dbl), allocatable           :: tmp1(:), tmp2(:)
    
    gin(1:step,1:2,1:nfour) => grid(1:2*step*nfour)
    gout(1:step,1:nfour)    => grid(1:  step*nfour)
    
    allocate( tmp1(step), tmp2(step) )
    
    do i1 = 1, nfour
      tmp1 = gin(1:step,1,i1)
      tmp2 = gin(1:step,2,i1)
      
      do concurrent ( i2 = 1:step )
        gout(i2,i1) = tmp1(i2) * tmp2(i2)
      end do
    end do
    
    deallocate( tmp1, tmp2 )
    
  end subroutine grid_op_vcss_sub
  
  module pure subroutine vcss_sub(this, cajm, cbjm, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: cc(:), cr(:)
    
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
    
  end subroutine vcss_sub
  
end submodule vcss